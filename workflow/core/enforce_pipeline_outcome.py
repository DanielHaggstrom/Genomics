#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys

from workflow.core.pipeline_lib import load_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    summary = load_json(args.summary)
    failed = summary.get("failed_enabled_stages") or []
    if failed:
        print(
            "Pipeline completed with stage failures: "
            + ", ".join(failed)
            + f". See {summary['report_path']} for diagnostics.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    print(f"Pipeline completed successfully. Report: {summary['report_path']}")


if __name__ == "__main__":
    main()
