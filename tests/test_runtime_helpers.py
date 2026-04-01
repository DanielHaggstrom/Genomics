from __future__ import annotations

import unittest

from workflow.core.runtime import (
    format_java_memory_mb,
    parse_memory_to_gb,
    parse_memory_to_mb,
    plan_contig_shards,
    plan_worker_pool,
)


class RuntimeHelperTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
