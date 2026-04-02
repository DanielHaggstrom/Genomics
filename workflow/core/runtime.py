#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    with open(destination, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def ensure_dir(path: str | pathlib.Path) -> pathlib.Path:
    destination = pathlib.Path(path)
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def read_text_if_exists(path: str | pathlib.Path) -> str | None:
    target = pathlib.Path(path)
    if not target.exists():
        return None
    return target.read_text(encoding="utf-8")


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


def artifacts_exist(paths: Iterable[str | None]) -> bool:
    for raw in paths:
        if not raw:
            continue
        if not pathlib.Path(raw).exists():
            return False
    return True


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


def complete_artifact(main_path: str | pathlib.Path, *, extra_sidecars: Iterable[str | pathlib.Path] | None = None) -> bool:
    target = pathlib.Path(main_path)
    if not target.exists():
        return False
    sidecars = []
    inferred = infer_index_path(target)
    if inferred is not None:
        sidecars.append(inferred)
    if extra_sidecars:
        sidecars.extend(pathlib.Path(p) for p in extra_sidecars)
    return all(sidecar.exists() for sidecar in sidecars)


def run_command(
    args: list[str],
    *,
    cwd: str | pathlib.Path | None = None,
    env: dict[str, str] | None = None,
    log_path: str | pathlib.Path | None = None,
    check: bool = True,
    timeout: int | float | None = None,
) -> subprocess.CompletedProcess[str]:
    final_env = os.environ.copy()
    if env:
        final_env.update(env)

    completed = subprocess.run(
        args,
        cwd=str(cwd) if cwd else None,
        env=final_env,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )

    if log_path is not None:
        log = pathlib.Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "w", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(args) + "\n")
            if completed.stdout:
                handle.write(completed.stdout)
                if not completed.stdout.endswith("\n"):
                    handle.write("\n")
            if completed.stderr:
                handle.write("\n[stderr]\n")
                handle.write(completed.stderr)
                if not completed.stderr.endswith("\n"):
                    handle.write("\n")

    if check and completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode,
            args,
            output=completed.stdout,
            stderr=completed.stderr,
        )

    return completed


def remove_if_exists(*paths: str | pathlib.Path) -> None:
    for raw in paths:
        target = pathlib.Path(raw)
        if target.exists() or target.is_symlink():
            target.unlink()


def atomic_copy(source: str | pathlib.Path, destination: str | pathlib.Path) -> pathlib.Path:
    src = pathlib.Path(source)
    dst = pathlib.Path(destination)
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
    destination_dir = ensure_dir(destination_dir)
    destination = destination_dir / source.name
    metadata_path = destination_dir / f"{source.name}.staged.json"
    fingerprint = fingerprint or sample_file_fingerprint(source)
    if destination.exists() and metadata_path.exists():
        try:
            metadata = load_json(metadata_path)
        except Exception:
            metadata = {}
        if metadata.get("source_fingerprint", {}).get("id") == fingerprint["id"]:
            return str(destination)

    remove_if_exists(destination)
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
    mounts: dict[str | pathlib.Path, str],
    inner_args: list[str],
    docker_bin: str,
    env: dict[str, str] | None = None,
) -> list[str]:
    command = [docker_bin, "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}"]
    for host_path, container_path in mounts.items():
        command.extend(["-v", f"{docker_mount_path(host_path, docker_bin)}:{container_path}"])
    if env:
        for key, value in env.items():
            command.extend(["-e", f"{key}={value}"])
    command.append(image)
    command.extend(inner_args)
    return command


def format_java_memory(raw: str) -> str:
    lowered = raw.lower().replace(" ", "")
    lowered = lowered.replace(".gb", "g").replace("gb", "g")
    lowered = lowered.replace(".mb", "m").replace("mb", "m")
    return lowered


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
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {executor.submit(task): idx for idx, task in enumerate(task_list)}
        for future in as_completed(future_map):
            idx = future_map[future]
            result = future.result()
            results[idx] = result
            if on_complete is not None:
                on_complete(idx, result)
    return results
