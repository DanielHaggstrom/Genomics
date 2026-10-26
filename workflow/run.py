"""Local stage driver. Successful checkpoints are reused on every invocation."""
from __future__ import annotations

import argparse
import fcntl
import pathlib
import sys

from workflow.core.manifest import cmd_validate, resolve_linux_cache_root
from workflow.core.reference import cmd_prepare_context, project_cache_id, validate_work_roots
from workflow.core.runtime import REPO_ROOT, interruptible_commands, load_json, write_json
from workflow.reporting.report import enforce_outcome, render_report
from workflow.stages import monogenic, prs


def run(manifest_path: pathlib.Path) -> int:
    manifest_path = manifest_path.resolve()
    cache_root = pathlib.Path(resolve_linux_cache_root(manifest_path))
    run_dir = cache_root / project_cache_id() / "run"
    # Check the actual destination too, in case an existing cache child is a symlink.
    from workflow.core.runtime import load_yaml, normalize_path
    raw = load_yaml(manifest_path)
    validate_work_roots(
        runtime_root=run_dir, results_root=run_dir,
        input_paths=[normalize_path(raw["inputs"][key], base_dir=manifest_path.parent) for key in ("cram", "crai")],
        protected_roots=[normalize_path(path, base_dir=manifest_path.parent) for path in raw["inputs"].get("protected_roots", [])],
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir.parent / "run.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("Another analysis is already using this project cache.", file=sys.stderr)
            return 1
        phase = "validation"

        def status(state: str, message: str | None = None) -> None:
            write_json(run_dir / "status.json", {
                "status": state, "phase": phase, "manifest_path": str(manifest_path), "message": message,
            })

        try:
            with interruptible_commands():
                status("running")
                validated = run_dir / "validated_manifest.json"
                cmd_validate(argparse.Namespace(manifest=str(manifest_path),
                             schema=str(REPO_ROOT / "manifests/schema/run_manifest.schema.json"), output=str(validated)))
                phase = "preparation"
                status("running")
                context_path = run_dir / "context.json"
                cmd_prepare_context(argparse.Namespace(validated_manifest=str(validated), output=str(context_path)))
                context = load_json(context_path)
                stages = {}
                for name, module in (("monogenic", monogenic), ("prs", prs)):
                    destination = run_dir / f"{name}_stage.json"
                    if context["manifest"]["analysis"][name]["enabled"]:
                        phase = name
                        status("running")
                        print(f"Stage: {name}", flush=True)
                        module.main(["--context", str(context_path), "--output", str(destination)])
                        stages[name] = load_json(destination)
                    else:
                        stages[name] = {"stage": name, "status": "skipped"}
                        write_json(destination, stages[name])
                phase = "report"
                status("running")
                summary = run_dir / "summary.json"
                render_report(context, stages["monogenic"], stages["prs"], output_summary=summary)
                result = load_json(summary)
                enforce_outcome(result)
                status("ok")
                print(f"Report: {result['report_path']}")
                return 0
        except KeyboardInterrupt:
            status("stopped", "Rerun the same command to reuse completed checkpoints.")
            print("Stopped. Completed checkpoints retained; rerun the same command to resume.", file=sys.stderr)
            return 130
        except (Exception, SystemExit) as exc:
            status("failed", str(exc))
            print(f"Failed during {phase}: {exc}", file=sys.stderr)
            return 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument("--resume", "-resume", action="store_true", help="Compatibility option; resume is always enabled")
    args = parser.parse_args()
    try:
        return run(args.manifest)
    except (Exception, SystemExit) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
