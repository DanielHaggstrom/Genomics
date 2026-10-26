from __future__ import annotations

import unittest
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from workflow.core import runtime

from workflow.core.runtime import (
    atomic_copy,
    build_docker_run_command,
    complete_artifact,
    complete_artifacts_exist,
    checkpoint_complete,
    format_java_memory_mb,
    init_shard_progress,
    germline_gvcf_dir,
    interruptible_commands,
    parse_memory_to_gb,
    parse_memory_to_mb,
    plan_contig_shards,
    plan_worker_pool,
    record_shard_completion,
    record_checkpoint,
    run_command,
    run_tasks_in_parallel,
    stage_file,
    write_json,
)


class RuntimeHelperTests(unittest.TestCase):
    def test_cancelled_peer_cannot_hide_the_original_failure(self) -> None:
        def failure():
            raise RuntimeError("original shard error")

        # Deliver the cancelled peer first, regardless of thread scheduling.
        with patch.object(runtime, "as_completed", side_effect=lambda futures: reversed(list(futures))), \
             self.assertRaisesRegex(RuntimeError, "original shard error"), interruptible_commands():
            run_tasks_in_parallel([failure, runtime.check_cancelled], max_workers=2)

    def test_download_resumes_partial_and_publishes_only_valid_success(self) -> None:
        with TemporaryDirectory() as directory:
            destination = Path(directory) / "download.gz"
            partial = Path(str(destination) + ".part")

            def interrupted(args, **kwargs):
                partial.write_bytes(b"partial")
                raise KeyboardInterrupt

            with patch.object(runtime, "run_command", side_effect=interrupted), self.assertRaises(KeyboardInterrupt):
                runtime.download_url_to_path("https://example.invalid/fixture", destination)
            self.assertFalse(destination.exists())

            def resumed(args, **kwargs):
                self.assertEqual(args[args.index("--continue-at") + 1], "-")
                self.assertEqual(partial.read_bytes(), b"partial")
                partial.write_bytes(b"complete")

            with patch.object(runtime, "run_command", side_effect=resumed):
                runtime.download_url_to_path("https://example.invalid/fixture", destination,
                                             validator=lambda path: path.read_bytes() == b"complete")
            self.assertEqual(destination.read_bytes(), b"complete")
            self.assertFalse(partial.exists())
            with patch.object(runtime, "run_command", side_effect=lambda *a, **k: partial.write_bytes(b"invalid")), \
                 self.assertRaisesRegex(RuntimeError, "failed validation"):
                runtime.download_url_to_path("https://example.invalid/fixture", destination, validator=lambda path: False)
            self.assertEqual(destination.read_bytes(), b"complete")
            self.assertFalse(partial.exists())

    def test_json_update_preserves_last_valid_state_on_failed_write(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "state.json"
            write_json(path, {"completed": ["chr1"]})
            with self.assertRaises(TypeError):
                write_json(path, {"invalid": object()})
            self.assertEqual(json.loads(path.read_text()), {"completed": ["chr1"]})
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_checkpoint_requires_success_receipt_and_unchanged_outputs(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data, index, receipt = root / "data.vcf.gz", root / "data.vcf.gz.tbi", root / "completed.json"
            data.write_bytes(b"completed fixture")
            index.write_bytes(b"index")
            self.assertFalse(checkpoint_complete(receipt, [data, index]))
            record_checkpoint(receipt, [data, index])
            self.assertTrue(checkpoint_complete(receipt, [data, index]))
            data.write_bytes(b"interrupted replacement")
            self.assertFalse(checkpoint_complete(receipt, [data, index]))
            receipt.write_text("{partial")
            self.assertFalse(checkpoint_complete(receipt, [data, index]))
            with self.assertRaises(ValueError):
                record_checkpoint(receipt, [root / "missing"])

    def test_logged_commands_stream_output_before_exit_and_bound_capture(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ready, release, log = root / "ready", root / "release", root / "command.log"
            script = ("import pathlib,sys,time; print('live output',flush=True); "
                      "pathlib.Path(sys.argv[1]).touch(); "
                      "exec(\"while not pathlib.Path(sys.argv[2]).exists(): time.sleep(0.01)\"); "
                      "print('x'*200000)")
            results = []
            worker = threading.Thread(target=lambda: results.append(run_command([sys.executable, "-c", script, str(ready), str(release)], log_path=log)))
            worker.start()
            try:
                deadline = time.monotonic() + 3
                while not ready.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertTrue(ready.exists())
                self.assertIn("live output", log.read_text())
                self.assertTrue(worker.is_alive())
            finally:
                release.touch()
                worker.join(timeout=3)
            self.assertFalse(worker.is_alive())
            self.assertEqual(results[0].stdout, "")
            self.assertLessEqual(len(results[0].stderr), 65536)
            self.assertGreater(log.stat().st_size, 200000)
            captured = run_command([sys.executable, "-c", "print('header')"], log_path=log, capture_output=True)
            self.assertEqual(captured.stdout, "header\n")

    def test_failed_shard_stops_active_process_group_and_skips_queued_work(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pid_file, queued = root / "pid", root / "queued"
            script = ("import pathlib,subprocess,sys,time; "
                      "child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)']); "
                      "pathlib.Path(sys.argv[1]).write_text(str(child.pid)); time.sleep(30)")
            def fail_after_start():
                deadline = time.monotonic() + 3
                while not pid_file.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                if not pid_file.exists():
                    raise AssertionError("Synthetic child did not start")
                raise RuntimeError("synthetic shard failed")
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, "synthetic shard failed"), interruptible_commands():
                run_tasks_in_parallel([
                    lambda: run_command([sys.executable, "-c", script, str(pid_file)], log_path=root / "child.log"),
                    fail_after_start,
                    lambda: queued.touch(),
                ], max_workers=2)
            self.assertLess(time.monotonic() - started, 5)
            self.assertFalse(queued.exists())
            child_stat = Path(f"/proc/{pid_file.read_text()}/stat")
            self.assertTrue(not child_stat.exists() or child_stat.read_text().split()[2] == "Z")

    def test_sigterm_stops_command_and_exits_without_leaving_child_running(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pid_file = root / "pid"
            command = ("import pathlib,os,sys,time; pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); "
                       "time.sleep(30)")
            driver = ("import sys; from workflow.core.runtime import interruptible_commands,run_command\n"
                      "try:\n"
                      " with interruptible_commands():\n"
                      f"  run_command([sys.executable, '-c', {command!r}, sys.argv[1]])\n"
                      "except KeyboardInterrupt: sys.exit(130)\n")
            process = subprocess.Popen([sys.executable, "-c", driver, str(pid_file)], start_new_session=True,
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            try:
                deadline = time.monotonic() + 3
                while not pid_file.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertTrue(pid_file.exists())
                process.send_signal(signal.SIGTERM)
                self.assertEqual(process.wait(timeout=5), 130)
                with self.assertRaises(ProcessLookupError):
                    os.kill(int(pid_file.read_text()), 0)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()

    def test_required_artifacts_reject_absent_empty_or_directory_values(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            empty = tmp / "empty.txt"
            empty.touch()
            for paths in ([], [None], [""], [tmp], [empty]):
                with self.subTest(paths=paths):
                    self.assertFalse(complete_artifacts_exist(paths))
            vcf = tmp / "sample.vcf.gz"
            vcf.write_bytes(b"synthetic data")
            index = Path(str(vcf) + ".tbi")
            index.touch()
            self.assertFalse(complete_artifact(vcf))

    def test_copy_and_staging_refuse_source_aliases_without_changing_them(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            original = tmp / "original.txt"
            original.write_bytes(b"preserve source")
            alias = tmp / "alias.txt"
            alias.hardlink_to(original)
            with self.assertRaisesRegex(ValueError, "source"):
                atomic_copy(original, alias)
            with self.assertRaisesRegex(ValueError, "source"):
                stage_file(original, tmp)
            self.assertEqual(original.read_bytes(), b"preserve source")
            self.assertEqual(alias.read_bytes(), b"preserve source")

    def test_docker_mount_modes_keep_input_aliases_read_only(self) -> None:
        command = build_docker_run_command(
            image="test-image", docker_bin="docker", inner_args=["true"],
            mounts=[("/data", "/input", "ro"), ("/data", "/ref", "ro"), ("/cache", "/out", "rw")],
        )
        self.assertIn("/data:/input:ro", command)
        self.assertIn("/data:/ref:ro", command)
        self.assertIn("/cache:/out:rw", command)
        with self.assertRaisesRegex(ValueError, "mount mode"):
            build_docker_run_command(image="test", docker_bin="docker", inner_args=[], mounts={"/data": ("/input", "invalid")})

    def test_shared_gvcf_cache_changes_with_reference_or_calling_inputs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cram, reference = tmp / "synthetic.cram", tmp / "synthetic.fa"
            cram.write_bytes(b"synthetic alignment")
            reference.write_bytes(b">chr1\nACGT\n")
            args = dict(cram_path=cram, reference_fasta=reference, image="gatk:test", intervals=["chr1"])
            first = germline_gvcf_dir(tmp / "cache", **args)
            self.assertEqual(first, germline_gvcf_dir(tmp / "cache", **args))
            self.assertNotEqual(first, germline_gvcf_dir(tmp / "cache", **{**args, "intervals": ["chr2"]}))
            reference.write_bytes(b">chr1\nACGA\n")
            self.assertNotEqual(first, germline_gvcf_dir(tmp / "cache", **args))

    def test_complete_artifacts_exist_requires_vcf_sidecars(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf = tmp / "sample.vcf.gz"
            vcf.write_text("placeholder", encoding="utf-8")
            self.assertFalse(complete_artifacts_exist([vcf]))
            (tmp / "sample.vcf.gz.tbi").write_text("index", encoding="utf-8")
            self.assertTrue(complete_artifacts_exist([vcf]))

    def test_plan_contig_shards_groups_non_primary_contigs_into_alt(self) -> None:
        shards = plan_contig_shards(["chr1", "chr2", "chrUn_KI270442v1", "chrEBV"], "chr")
        self.assertEqual(shards[0]["name"], "chr1")
        self.assertEqual(shards[1]["name"], "chr2")
        self.assertEqual(shards[2]["name"], "alt")
        self.assertEqual(shards[2]["intervals"], ["chrUn_KI270442v1", "chrEBV"])

    def test_memory_helpers_accept_manifest_style_values(self) -> None:
        self.assertEqual(parse_memory_to_mb("14.GB"), 14 * 1024)
        self.assertAlmostEqual(parse_memory_to_gb("1536mb"), 1.5)
        self.assertEqual(format_java_memory_mb(2048), "2048m")

    def test_plan_worker_pool_scales_automatically(self) -> None:
        plan = plan_worker_pool(
            stage_name="monogenic",
            shard_count=26,
            max_cpus=8,
            max_memory="14.GB",
            default_worker_cap=4,
            default_pairhmm_cap=2,
            min_memory_per_worker_gb=3.0,
            free_gb=100.0,
            min_free_gb=40.0,
        )
        self.assertEqual(plan["workers"], 4)
        self.assertEqual(plan["pairhmm_threads"], 2)
        self.assertEqual(plan["java_memory"], "3328m")
        self.assertEqual(plan["worker_source"], "automatic")
        self.assertEqual(plan["pairhmm_source"], "automatic")

    def test_plan_worker_pool_caps_requested_overrides(self) -> None:
        plan = plan_worker_pool(
            stage_name="prs_dense_genotypes",
            shard_count=26,
            max_cpus=8,
            max_memory="14.GB",
            requested_workers=10,
            requested_pairhmm_threads=8,
            default_worker_cap=4,
            default_pairhmm_cap=2,
            min_memory_per_worker_gb=2.5,
            free_gb=100.0,
            min_free_gb=40.0,
        )
        self.assertEqual(plan["workers"], 5)
        self.assertEqual(plan["pairhmm_threads"], 1)
        self.assertEqual(plan["worker_source"], "override")
        self.assertEqual(plan["pairhmm_source"], "override")
        self.assertIn("worker_override_note", plan)
        self.assertIn("pairhmm_override_note", plan)

    def test_plan_worker_pool_reduces_workers_when_disk_headroom_is_low(self) -> None:
        plan = plan_worker_pool(
            stage_name="monogenic",
            shard_count=26,
            max_cpus=8,
            max_memory="14.GB",
            default_worker_cap=4,
            default_pairhmm_cap=2,
            min_memory_per_worker_gb=3.0,
            free_gb=45.0,
            min_free_gb=40.0,
        )
        self.assertEqual(plan["workers"], 2)
        self.assertEqual(plan["pairhmm_threads"], 2)
        self.assertEqual(plan["limits"]["disk_worker_limit"], 2)

    def test_run_tasks_in_parallel_reports_completed_tasks(self) -> None:
        seen: list[tuple[int, str]] = []
        results = run_tasks_in_parallel(
            [lambda: "first", lambda: "second"],
            max_workers=2,
            on_complete=lambda idx, result: seen.append((idx, result)),
        )
        self.assertEqual(results, ["first", "second"])
        self.assertEqual(sorted(seen), [(0, "first"), (1, "second")])

    def test_shard_progress_tracks_reused_and_rebuilt_counts(self) -> None:
        progress = init_shard_progress(phase="calling_shards", total_shards=3, score_hash="abc")
        record_shard_completion(progress, shard_name="chr1", reused=True)
        record_shard_completion(progress, shard_name="chr2", reused=False)
        self.assertEqual(progress["completed_shards"], 2)
        self.assertEqual(progress["pending_shards"], 1)
        self.assertEqual(progress["reused_completed_shards"], 1)
        self.assertEqual(progress["rebuilt_completed_shards"], 1)
        self.assertEqual(progress["last_completed_shard"], "chr2")


if __name__ == "__main__":
    unittest.main()
