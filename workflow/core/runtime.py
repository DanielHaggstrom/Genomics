#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
import uuid
from concurrent.futures import CancelledError, ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from typing import Any, Callable, Iterable

import yaml
from jsonschema import Draft202012Validator


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
WINDOWS_PATH_RE = re.compile(r"^[A-Za-z]:[\\/]")
MEMORY_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([gGmM][bB]?|\.[gG][bB]?|\.[mM][bB]?)\s*$")


def load_yaml(path: str | pathlib.Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_json(path: str | pathlib.Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: str | pathlib.Path, payload: Any) -> None:
    destination = pathlib.Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=destination.parent, delete=False) as handle:
            temporary = pathlib.Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def ensure_dir(path: str | pathlib.Path) -> pathlib.Path:
    destination = pathlib.Path(path)
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def is_windows_path(raw: str) -> bool:
    return bool(WINDOWS_PATH_RE.match(raw))


def is_mnt_path(raw: str | pathlib.Path) -> bool:
    return str(raw).startswith("/mnt/")


def normalize_path(raw: str | None, *, base_dir: pathlib.Path | None = None) -> str | None:
    if raw is None:
        return None

    raw = str(raw).strip()
    if not raw:
        return None

    if raw.startswith("~"):
        return str(pathlib.Path(raw).expanduser().resolve())

    if is_windows_path(raw):
        wslpath = shutil.which("wslpath")
        if wslpath:
            completed = subprocess.run(
                [wslpath, "-u", raw],
                check=True,
                capture_output=True,
                text=True,
            )
            return completed.stdout.strip()

        drive = raw[0].lower()
        rest = raw[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{drive}/{rest}"

    candidate = pathlib.Path(raw)
    if not candidate.is_absolute():
        candidate = (base_dir or REPO_ROOT).joinpath(candidate)
    return str(candidate.resolve())


def detect_executable(name: str) -> str | None:
    return shutil.which(name)


def validate_payload(payload: dict[str, Any], schema: dict[str, Any]) -> list[str]:
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(payload), key=lambda e: list(e.absolute_path))
    messages = []
    for error in errors:
        location = ".".join(str(x) for x in error.absolute_path) or "<root>"
        messages.append(f"{location}: {error.message}")
    return messages


def stage_payload(stage: str, status: str, **extra: Any) -> dict[str, Any]:
    payload = {"stage": stage, "status": status}
    payload.update(extra)
    return payload


def hash_strings(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def sample_file_fingerprint(path: str | pathlib.Path, *, chunk_size: int = 1024 * 1024) -> dict[str, Any]:
    target = pathlib.Path(path)
    stat = target.stat()
    digest = hashlib.sha256()
    digest.update(str(target).encode("utf-8"))
    digest.update(str(stat.st_size).encode("utf-8"))
    digest.update(str(stat.st_mtime_ns).encode("utf-8"))
    with open(target, "rb") as handle:
        first = handle.read(chunk_size)
        digest.update(first)
        if stat.st_size > chunk_size:
            handle.seek(max(stat.st_size - chunk_size, 0))
            digest.update(handle.read(chunk_size))
    return {
        "path": str(target),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "id": digest.hexdigest(),
        "method": "sha256(path,size,mtime,first_chunk,last_chunk)",
    }


def germline_gvcf_dir(
    cache_root: str | pathlib.Path, *, cram_path: str | pathlib.Path,
    reference_fasta: str | pathlib.Path, image: str, intervals: Iterable[str],
) -> pathlib.Path:
    """Shared whole-interval HaplotypeCaller output, independent of selected scores."""
    key = hash_strings([
        "germline-gvcf-v1", sample_file_fingerprint(cram_path)["id"],
        sample_file_fingerprint(reference_fasta)["id"], image, *intervals,
    ])
    return ensure_dir(pathlib.Path(cache_root) / key)


def artifacts_exist(paths: Iterable[str | pathlib.Path | None]) -> bool:
    required = list(paths)
    return bool(required) and all(
        raw and pathlib.Path(raw).is_file() and pathlib.Path(raw).stat().st_size > 0
        for raw in required
    )


def infer_index_path(path: str | pathlib.Path) -> pathlib.Path | None:
    target = pathlib.Path(path)
    name = target.name
    if name.endswith(".vcf.gz"):
        return pathlib.Path(str(target) + ".tbi")
    if name.endswith(".bcf"):
        return pathlib.Path(str(target) + ".csi")
    if name.endswith(".cram"):
        return pathlib.Path(str(target) + ".crai")
    return None


def complete_artifact(main_path: str | pathlib.Path | None, *, extra_sidecars: Iterable[str | pathlib.Path] | None = None) -> bool:
    if not main_path:
        return False
    target = pathlib.Path(main_path)
    if not target.is_file() or target.stat().st_size == 0:
        return False
    sidecars = []
    inferred = infer_index_path(target)
    if inferred is not None:
        sidecars.append(inferred)
    if extra_sidecars:
        sidecars.extend(pathlib.Path(p) for p in extra_sidecars)
    return all(sidecar.is_file() and sidecar.stat().st_size > 0 for sidecar in sidecars)


def complete_artifacts_exist(paths: Iterable[str | pathlib.Path | None]) -> bool:
    required = list(paths)
    return bool(required) and all(complete_artifact(raw) for raw in required)


def checkpoint_complete(receipt: str | pathlib.Path, artifacts: Iterable[str | pathlib.Path]) -> bool:
    paths = list(artifacts)
    try:
        return artifacts_exist(paths) and load_json(receipt) == {
            "artifacts": [sample_file_fingerprint(path) for path in paths],
        }
    except (OSError, ValueError, TypeError):
        return False


def record_checkpoint(receipt: str | pathlib.Path, artifacts: Iterable[str | pathlib.Path]) -> None:
    paths = list(artifacts)
    if not artifacts_exist(paths):
        raise ValueError("Cannot checkpoint missing or empty outputs")
    write_json(receipt, {"artifacts": [sample_file_fingerprint(path) for path in paths]})


_ACTIVE_COMMANDS: dict[subprocess.Popen, tuple[str, str] | None] = {}
_COMMAND_LOCK = threading.Lock()
_CANCEL_COMMANDS = threading.Event()


def check_cancelled() -> None:
    if _CANCEL_COMMANDS.is_set():
        raise CancelledError("Analysis cancelled")


def _stop_command(process: subprocess.Popen, container: tuple[str, str] | None) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    # Docker containers belong to the daemon, not the client's process group.
    # Each run has a unique name, so cancellation can remove only its own container.
    if container is not None:
        try:
            subprocess.run([container[0], "rm", "-f", container[1]], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=5, check=False)
        except (OSError, subprocess.TimeoutExpired):
            pass
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def cancel_commands() -> None:
    with _COMMAND_LOCK:
        _CANCEL_COMMANDS.set()
        commands = list(_ACTIVE_COMMANDS.items())
    for process, container in commands:
        _stop_command(process, container)


@contextmanager
def interruptible_commands():
    """One analysis run: turn stop signals into cleanup of its active commands."""
    def interrupted(_signum, _frame):
        raise KeyboardInterrupt("Analysis interrupted")
    previous = {signum: signal.getsignal(signum) for signum in (signal.SIGINT, signal.SIGTERM)}
    for signum in previous:
        signal.signal(signum, interrupted)
    _CANCEL_COMMANDS.clear()
    try:
        yield
    except BaseException:
        cancel_commands()
        raise
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
        _CANCEL_COMMANDS.clear()


def run_command(
    args: list[str],
    *,
    cwd: str | pathlib.Path | None = None,
    env: dict[str, str] | None = None,
    log_path: str | pathlib.Path | None = None,
    check: bool = True,
    timeout: int | float | None = None,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run a cancellable command; long-running command output goes straight to disk.

    Without a log path, output is captured as before. Logged commands return only
    a bounded diagnostic tail unless capture_output=True requests their full stdout.
    """
    final_env = {**os.environ, **(env or {})}
    command = list(args)
    container = None
    if pathlib.Path(command[0]).name in {"docker", "docker.exe"} and command[1:2] == ["run"]:
        container = (command[0], "genomics-" + uuid.uuid4().hex)
        command[2:2] = ["--name", container[1]]
    log = None
    process = None
    try:
        if log_path is not None:
            path = pathlib.Path(log_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            log = path.open("w", encoding="utf-8")
            log.write("$ " + " ".join(command) + "\n")
            log.flush()
        with _COMMAND_LOCK:
            check_cancelled()
            process = subprocess.Popen(
                command, cwd=str(cwd) if cwd else None, env=final_env, text=True,
                stdout=subprocess.PIPE if log is None or capture_output else log,
                stderr=subprocess.PIPE if log is None else log, start_new_session=True,
            )
            _ACTIVE_COMMANDS[process] = container
        stdout, stderr = process.communicate(timeout=timeout)
        if log is not None:
            if stdout:
                log.write(stdout)
            log.flush()
            with pathlib.Path(log_path).open("rb") as tail:
                tail.seek(max(0, tail.seek(0, os.SEEK_END) - 65536))
                stderr = tail.read().decode("utf-8", errors="replace")
        completed = subprocess.CompletedProcess(args, process.returncode, stdout or "", stderr or "")
        if check and completed.returncode != 0:
            raise subprocess.CalledProcessError(completed.returncode, args, output=completed.stdout, stderr=completed.stderr)
        return completed
    except BaseException:
        if process is not None:
            _stop_command(process, container)
        raise
    finally:
        if process is not None:
            with _COMMAND_LOCK:
                _ACTIVE_COMMANDS.pop(process, None)
        if log is not None:
            log.close()


def curl_args() -> list[str]:
    return ["curl", "-fsSL", "--retry", "4", "--retry-delay", "2",
            "--connect-timeout", "30", "--max-time", "21600"]


def fetch_url_text(url: str) -> str:
    return run_command([*curl_args(), url], check=True).stdout


def download_url_to_path(
    url: str, destination: pathlib.Path,
    *, validator: Callable[[pathlib.Path], bool] | None = None,
) -> pathlib.Path:
    ensure_dir(destination.parent)
    partial = destination.with_name(destination.name + ".part")
    # Curl resumes an interrupted transfer; only a successful, validated transfer
    # becomes the canonical cache file. Interruptions leave .part for the next run.
    run_command([*curl_args(), "--continue-at", "-", "-o", str(partial), url], check=True)
    if not partial.is_file() or partial.stat().st_size == 0 or (validator and not validator(partial)):
        remove_if_exists(partial)
        raise RuntimeError(f"Downloaded file failed validation: {url}")
    partial.replace(destination)
    return destination


def completed_vcf(path: pathlib.Path) -> bool:
    return checkpoint_complete(path.with_name(path.name + ".completed.json"),
                               [path, pathlib.Path(str(path) + ".tbi")])


def record_vcf_checkpoint(path: pathlib.Path) -> None:
    record_checkpoint(path.with_name(path.name + ".completed.json"),
                      [path, pathlib.Path(str(path) + ".tbi")])


def remove_if_exists(*paths: str | pathlib.Path) -> None:
    for raw in paths:
        target = pathlib.Path(raw)
        if target.exists() or target.is_symlink():
            target.unlink()


def atomic_copy(source: str | pathlib.Path, destination: str | pathlib.Path) -> pathlib.Path:
    src = pathlib.Path(source)
    dst = pathlib.Path(destination)
    if src.resolve() == dst.resolve() or (dst.exists() and src.samefile(dst)):
        raise ValueError(f"Refusing to overwrite source file while copying: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=str(dst.parent), delete=False) as handle:
        temp_path = pathlib.Path(handle.name)
    try:
        shutil.copy2(src, temp_path)
        temp_path.replace(dst)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return dst


def needs_staging(source_path: str | pathlib.Path, linux_cache_root: str | pathlib.Path) -> bool:
    source = pathlib.Path(source_path)
    linux_root = pathlib.Path(linux_cache_root)
    try:
        source.relative_to(linux_root)
        return False
    except ValueError:
        return is_mnt_path(source)


def stage_file(
    source_path: str | pathlib.Path,
    destination_dir: str | pathlib.Path,
    *,
    fingerprint: dict[str, Any] | None = None,
) -> str:
    source = pathlib.Path(source_path)
    destination_dir = pathlib.Path(destination_dir)
    destination = destination_dir / source.name
    if source.resolve() == destination.resolve() or (destination.exists() and source.samefile(destination)):
        raise ValueError(f"Refusing to stage a file over its source: {source}")
    ensure_dir(destination_dir)
    metadata_path = destination_dir / f"{source.name}.staged.json"
    fingerprint = fingerprint or sample_file_fingerprint(source)
    if destination.exists() and metadata_path.exists():
        try:
            metadata = load_json(metadata_path)
        except Exception:
            metadata = {}
        if metadata.get("source_fingerprint", {}).get("id") == fingerprint["id"]:
            return str(destination)

    atomic_copy(source, destination)
    write_json(
        metadata_path,
        {
            "source": str(source),
            "destination": str(destination),
            "source_fingerprint": fingerprint,
        },
    )
    return str(destination)


def primary_contigs(prefix: str) -> list[str]:
    chroms = [f"{prefix}{i}" for i in range(1, 23)]
    chroms.extend([f"{prefix}X", f"{prefix}Y"])
    chroms.extend([f"{prefix}M", f"{prefix}MT"])
    deduped = []
    seen = set()
    for chrom in chroms:
        if chrom in seen:
            continue
        seen.add(chrom)
        deduped.append(chrom)
    return deduped


def plan_contig_shards(contigs: Iterable[str], prefix: str) -> list[dict[str, Any]]:
    contig_list = list(contigs)
    primary = [name for name in primary_contigs(prefix) if name in contig_list]
    seen = set(primary)
    extras = [name for name in contig_list if name not in seen]
    shards = [{"name": name, "intervals": [name], "kind": "primary"} for name in primary]
    if extras:
        shards.append({"name": "alt", "intervals": extras, "kind": "alt"})
    return shards


def init_shard_progress(*, phase: str, total_shards: int, **extra: Any) -> dict[str, Any]:
    progress = {
        "phase": phase,
        "total_shards": int(total_shards),
        "completed_shards": 0,
        "pending_shards": int(total_shards),
        "reused_completed_shards": 0,
        "rebuilt_completed_shards": 0,
        "completed_shard_names": [],
    }
    progress.update(extra)
    return progress


def record_shard_completion(progress: dict[str, Any], *, shard_name: str, reused: bool) -> dict[str, Any]:
    completed = progress.setdefault("completed_shard_names", [])
    if shard_name not in completed:
        completed.append(shard_name)
        progress["completed_shards"] = len(completed)
        progress["pending_shards"] = max(0, int(progress.get("total_shards", len(completed))) - len(completed))
        if reused:
            progress["reused_completed_shards"] = int(progress.get("reused_completed_shards", 0)) + 1
        else:
            progress["rebuilt_completed_shards"] = int(progress.get("rebuilt_completed_shards", 0)) + 1
    progress["last_completed_shard"] = shard_name
    return progress


def ensure_docker_ready(docker_bin: str) -> None:
    completed = run_command([docker_bin, "version", "--format", "{{.Server.Version}}"], check=False)
    if completed.returncode == 0:
        return

    if docker_bin.endswith(".exe") and detect_executable("powershell.exe"):
        run_command(
            [
                "powershell.exe",
                "-NoProfile",
                "-Command",
                "& { $p='C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe'; if (Test-Path $p) { Start-Process -FilePath $p } }",
            ],
            check=False,
        )
        for _ in range(24):
            time.sleep(5)
            completed = run_command([docker_bin, "version", "--format", "{{.Server.Version}}"], check=False)
            if completed.returncode == 0:
                return

    raise SystemExit("Docker is not ready. Start Docker Desktop and retry.")


def docker_mount_path(path: str | pathlib.Path, docker_bin: str) -> str:
    raw = str(path)
    if docker_bin.endswith(".exe"):
        wslpath = detect_executable("wslpath")
        if wslpath:
            completed = subprocess.run([wslpath, "-w", raw], check=True, capture_output=True, text=True)
            return completed.stdout.strip()
    return raw


def build_docker_run_command(
    *,
    image: str,
    mounts: dict[str | pathlib.Path, str | tuple[str, str]] | list[tuple[str | pathlib.Path, str, str]],
    inner_args: list[str],
    docker_bin: str,
    env: dict[str, str] | None = None,
) -> list[str]:
    command = [docker_bin, "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}"]
    if isinstance(mounts, dict):
        entries = [
            (host, *target) if isinstance(target, tuple) else (host, target, "rw")
            for host, target in mounts.items()
        ]
    else:
        entries = mounts
    for host_path, container_path, mode in entries:
        if mode not in {"ro", "rw"}:
            raise ValueError(f"Unsupported Docker mount mode: {mode}")
        command.extend(["-v", f"{docker_mount_path(host_path, docker_bin)}:{container_path}:{mode}"])
    if env:
        for key, value in env.items():
            command.extend(["-e", f"{key}={value}"])
    command.append(image)
    command.extend(inner_args)
    return command


def parse_memory_to_mb(raw: str | int | float) -> int:
    if isinstance(raw, (int, float)):
        return int(raw)

    value = str(raw).strip()
    match = MEMORY_RE.match(value)
    if not match:
        raise ValueError(f"Unsupported memory value: {raw}")

    amount = float(match.group(1))
    unit = match.group(2).lower().replace(".", "")
    if unit in {"g", "gb"}:
        return int(amount * 1024)
    if unit in {"m", "mb"}:
        return int(amount)
    raise ValueError(f"Unsupported memory unit in value: {raw}")


def parse_memory_to_gb(raw: str | int | float) -> float:
    return parse_memory_to_mb(raw) / 1024.0


def format_java_memory_mb(memory_mb: int) -> str:
    return f"{int(memory_mb)}m"


def plan_worker_pool(
    *,
    stage_name: str,
    shard_count: int,
    max_cpus: int,
    max_memory: str,
    requested_workers: int | None = None,
    requested_pairhmm_threads: int | None = None,
    default_worker_cap: int,
    default_pairhmm_cap: int,
    min_memory_per_worker_gb: float,
    reserve_memory_gb: float = 1.0,
    target_cpus_per_worker: int = 2,
    free_gb: float | None = None,
    min_free_gb: float | None = None,
    low_disk_headroom_gb: float = 10.0,
) -> dict[str, Any]:
    shard_count = max(1, int(shard_count))
    max_cpus = max(1, int(max_cpus))
    total_memory_mb = max(1024, parse_memory_to_mb(max_memory))
    usable_memory_mb = max(1024, total_memory_mb - int(reserve_memory_gb * 1024))
    min_memory_per_worker_mb = max(1024, int(min_memory_per_worker_gb * 1024))

    cpu_limit_for_target = max(1, max_cpus // max(1, target_cpus_per_worker))
    memory_worker_limit = max(1, usable_memory_mb // min_memory_per_worker_mb)
    disk_worker_limit = shard_count
    disk_headroom_gb = None
    if free_gb is not None and min_free_gb is not None:
        disk_headroom_gb = max(0.0, float(free_gb) - float(min_free_gb))
        if disk_headroom_gb < low_disk_headroom_gb:
            disk_worker_limit = min(disk_worker_limit, 2)

    default_workers = min(shard_count, default_worker_cap, cpu_limit_for_target, memory_worker_limit, disk_worker_limit)
    hard_worker_limit = min(shard_count, max_cpus, memory_worker_limit, disk_worker_limit)
    if requested_workers is None:
        workers = max(1, default_workers)
        worker_source = "automatic"
    else:
        workers = max(1, min(int(requested_workers), hard_worker_limit))
        worker_source = "override"

    cpu_budget_per_worker = max(1, max_cpus // workers)
    default_pairhmm_threads = max(1, min(default_pairhmm_cap, cpu_budget_per_worker))
    hard_pairhmm_limit = max(1, min(cpu_budget_per_worker, default_pairhmm_cap))
    if requested_pairhmm_threads is None:
        pairhmm_threads = default_pairhmm_threads
        pairhmm_source = "automatic"
    else:
        pairhmm_threads = max(1, min(int(requested_pairhmm_threads), hard_pairhmm_limit))
        pairhmm_source = "override"

    java_memory_mb = max(1024, usable_memory_mb // workers)
    plan = {
        "stage_name": stage_name,
        "shard_count": shard_count,
        "workers": workers,
        "pairhmm_threads": pairhmm_threads,
        "java_memory_mb": java_memory_mb,
        "java_memory": format_java_memory_mb(java_memory_mb),
        "worker_source": worker_source,
        "pairhmm_source": pairhmm_source,
        "requested_workers": requested_workers,
        "requested_pairhmm_threads": requested_pairhmm_threads,
        "limits": {
            "default_worker_cap": default_worker_cap,
            "cpu_limit_for_target": cpu_limit_for_target,
            "memory_worker_limit": memory_worker_limit,
            "disk_worker_limit": disk_worker_limit,
            "hard_worker_limit": hard_worker_limit,
            "hard_pairhmm_limit": hard_pairhmm_limit,
        },
        "resources": {
            "max_cpus": max_cpus,
            "max_memory_mb": total_memory_mb,
            "usable_memory_mb": usable_memory_mb,
            "min_memory_per_worker_mb": min_memory_per_worker_mb,
            "free_gb": free_gb,
            "min_free_gb": min_free_gb,
            "disk_headroom_gb": disk_headroom_gb,
        },
    }
    if requested_workers is not None and int(requested_workers) != workers:
        plan["worker_override_note"] = f"Requested workers capped to {workers} by available resources."
    if requested_pairhmm_threads is not None and int(requested_pairhmm_threads) != pairhmm_threads:
        plan["pairhmm_override_note"] = f"Requested pairhmm_threads capped to {pairhmm_threads} by per-worker CPU budget."
    return plan


def run_tasks_in_parallel(
    tasks: Iterable[Callable[[], Any]],
    *,
    max_workers: int,
    on_complete: Callable[[int, Any], None] | None = None,
) -> list[Any]:
    task_list = list(tasks)
    if not task_list:
        return []
    if max_workers <= 1:
        results = []
        for idx, task in enumerate(task_list):
            result = task()
            results.append(result)
            if on_complete is not None:
                on_complete(idx, result)
        return results

    results: list[Any] = [None] * len(task_list)
    executor = ThreadPoolExecutor(max_workers=max_workers)
    failures: list[BaseException] = []
    failure_lock = threading.Lock()

    def execute(task: Callable[[], Any]) -> Any:
        check_cancelled()
        try:
            return task()
        except BaseException as exc:
            with failure_lock:
                if not failures and not isinstance(exc, CancelledError):
                    failures.append(exc)
                _CANCEL_COMMANDS.set()
            raise

    future_map = {}
    try:
        for idx, task in enumerate(task_list):
            future_map[executor.submit(execute, task)] = idx
        for future in as_completed(future_map):
            idx = future_map[future]
            result = future.result()
            results[idx] = result
            if on_complete is not None:
                on_complete(idx, result)
    except BaseException as exc:
        for future in future_map:
            future.cancel()
        cancel_commands()
        # A cancelled peer can finish before the shard that caused cancellation.
        if failures and not isinstance(exc, KeyboardInterrupt):
            raise failures[0]
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
        _CANCEL_COMMANDS.clear()
    return results
