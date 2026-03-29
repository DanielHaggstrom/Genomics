#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
from typing import Any

import yaml
from jsonschema import Draft202012Validator


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_yaml(path: str | pathlib.Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_json(path: str | pathlib.Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: str | pathlib.Path, payload: Any) -> None:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def ensure_dir(path: str | pathlib.Path) -> pathlib.Path:
    out = pathlib.Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def repo_path(*parts: str) -> pathlib.Path:
    return REPO_ROOT.joinpath(*parts)


WINDOWS_PATH_RE = re.compile(r"^[A-Za-z]:[\\/]")


def is_windows_path(raw: str) -> bool:
    return bool(WINDOWS_PATH_RE.match(raw))


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


def normalize_optional_paths(payload: dict[str, Any], *, base_dir: pathlib.Path | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in payload.items():
        if isinstance(value, dict):
            out[key] = normalize_optional_paths(value, base_dir=base_dir)
        elif isinstance(value, list):
            out[key] = value[:]
        elif isinstance(value, str) and ("/" in value or "\\" in value or value.startswith(".")):
            out[key] = normalize_path(value, base_dir=base_dir)
        else:
            out[key] = value
    return out


def run_command(
    args: list[str],
    *,
    cwd: str | pathlib.Path | None = None,
    env: dict[str, str] | None = None,
    log_path: str | pathlib.Path | None = None,
    check: bool = True,
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
    )

    if log_path is not None:
        log_path = pathlib.Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w", encoding="utf-8") as handle:
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


def sample_file_fingerprint(path: str | pathlib.Path, *, chunk_size: int = 1024 * 1024) -> dict[str, Any]:
    path = pathlib.Path(path)
    stat = path.stat()
    digest = hashlib.sha256()
    digest.update(str(path).encode("utf-8"))
    digest.update(str(stat.st_size).encode("utf-8"))
    digest.update(str(stat.st_mtime_ns).encode("utf-8"))
    with open(path, "rb") as handle:
        first = handle.read(chunk_size)
        digest.update(first)
        if stat.st_size > chunk_size:
            handle.seek(max(stat.st_size - chunk_size, 0))
            digest.update(handle.read(chunk_size))
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "id": digest.hexdigest(),
        "method": "sha256(path,size,mtime,first_chunk,last_chunk)",
    }


def hash_strings(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def validate_payload(payload: dict[str, Any], schema: dict[str, Any]) -> list[str]:
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(payload), key=lambda e: list(e.absolute_path))
    messages = []
    for error in errors:
        location = ".".join(str(x) for x in error.absolute_path) or "<root>"
        messages.append(f"{location}: {error.message}")
    return messages


def detect_executable(name: str) -> str | None:
    return shutil.which(name)


def stage_payload(stage: str, status: str, **extra: Any) -> dict[str, Any]:
    payload = {
        "stage": stage,
        "status": status,
    }
    payload.update(extra)
    return payload


def artifacts_exist(paths: list[str | None]) -> bool:
    for path in paths:
        if not path:
            continue
        if not pathlib.Path(path).exists():
            return False
    return True


def read_text_if_exists(path: str | pathlib.Path) -> str | None:
    path = pathlib.Path(path)
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8")
