#!/usr/bin/env python3
"""Content-aware QC checks used before final report generation."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from pipeline_common import analysis_run_dir, dock_run_dir, read_tsv, update_stage_status, write_tsv


def comparison_qc_rows(analysis: Path):
    root = analysis / "comparisons"
    groups = read_tsv(root / "receptor_ligand_groups.tsv")
    eligibility = read_tsv(root / "comparison_eligibility.tsv")
    results = read_tsv(root / "binding_mode_comparison.tsv")
    skips = read_tsv(root / "comparison_skips.tsv")
    single = sum(row.get("status") == "SKIPPED_SINGLE_LIGAND" for row in groups)
    multi = sum(int(row.get("unique_cas_count") or 0) >= 2 for row in groups)
    eligible_rows = [row for row in eligibility if row.get("eligible") == "1"]
    eligible_pairs = {
        (row.get("canonical_receptor_key"), row.get("cas_A"), row.get("cas_B"))
        for row in eligible_rows
    }
    successful = sum(row.get("comparison_status") == "SUCCESS" for row in results)
    partial = sum(row.get("comparison_status", "").startswith("PARTIAL") for row in results)
    result_ids = {row.get("comparison_id") for row in results if row.get("comparison_id")}
    runtime_skip_keys = {
        (row.get("canonical_receptor_key"), row.get("cas_A"), row.get("cas_B"), row.get("pocket_index"))
        for row in skips if row.get("status") in {"SKIPPED_CLUSTER_RESULT_MISSING", "FAILED_COMPARISON_EXCEPTION"}
    }
    accounted = set(result_ids)
    for row in eligible_rows:
        key = (row.get("canonical_receptor_key"), row.get("cas_A"), row.get("cas_B"), row.get("pocket_index"))
        if key in runtime_skip_keys:
            accounted.add(row.get("comparison_id"))
    unaccounted = {row.get("comparison_id") for row in eligible_rows} - accounted
    exceptions = sum(row.get("status") == "FAILED_COMPARISON_EXCEPTION" for row in skips)
    structural_skips = sum(row.get("status") in {
        "SKIPPED_RECEPTOR_STRUCTURE_MISMATCH", "SKIPPED_POCKET_DEFINITION_MISMATCH",
    } for row in skips)
    runtime_skips = sum(row.get("status") == "SKIPPED_CLUSTER_RESULT_MISSING" for row in skips)
    comparison_status = "FAIL" if exceptions or unaccounted else (
        "WARNING" if partial or runtime_skips else "PASS")
    skipped_status = "WARNING" if structural_skips or runtime_skips else "PASS"

    def row(check, observed, expected, status, message):
        return {"check": check, "status": status, "observed": observed, "expected": expected, "message": message}

    return [
        row("receptor_groups_total", len(groups), ">=0", "PASS", "casefolded receptor groups discovered"),
        row("single_ligand_receptor_groups", single, ">=0", "PASS", "single-CAS groups are valid and not comparable"),
        row("multi_ligand_receptor_groups", multi, ">=0", "PASS", "zero multi-ligand groups is not an error"),
        row("eligible_ligand_pairs", len(eligible_pairs), ">=0", "PASS", "unique eligible cross-CAS pairs"),
        row("eligible_pocket_comparisons", len(eligible_rows), ">=0", comparison_status,
            "eligible pocket comparisons; {} unaccounted and {} exceptions".format(len(unaccounted), exceptions)),
        row("successful_binding_mode_comparisons", successful, len(eligible_rows), comparison_status,
            "complete comparisons; partial and documented missing-cluster skips are accounted separately"),
        row("partial_binding_mode_comparisons", partial, 0, "WARNING" if partial else "PASS",
            "partial comparisons retain geometry/energy when interaction data are missing"),
        row("skipped_binding_mode_comparisons", len(skips), 0, skipped_status,
            "{} structural/grid skips; {} missing-cluster skips; single-ligand/no-common-pocket skips are valid".format(structural_skips, runtime_skips)),
    ]


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--analysis-run", type=Path, default=None); parser.add_argument("--dock-run", type=Path, default=None); args = parser.parse_args()
    analysis, dock = args.analysis_run or analysis_run_dir(), args.dock_run or dock_run_dir()
    qc = []
    poses = [r for r in read_tsv(analysis / "docking_post_analysis" / "pose_inventory.tsv") if r["extraction_status"] == "SUCCESS"]
    by_group = Counter((r["cas"], r["receptor"], r["pocket"]) for r in poses)
    for key in sorted({(r["cas"], r["receptor"], r["pocket"]) for r in read_tsv(dock / "job_manifest.tsv")}):
        count = by_group[key]
        qc.append({"check": "EXPECTED_100_POSES", "cas": key[0], "receptor": key[1], "pocket": key[2],
                   "status": "PASS" if count == 100 else "FAIL", "observed": count, "expected": 100,
                   "message": "validated extracted poses"})
    shard_results = read_tsv(dock / "docking_shard_results.tsv")
    failed_shards = sum(r["status"] not in {"SUCCESS", "SKIPPED_COMPLETE"} for r in shard_results)
    qc.append({"check": "ALL_DOCKING_SHARDS", "status": "PASS" if shard_results and failed_shards == 0 else "FAIL",
               "observed": failed_shards, "expected": 0, "message": "failed or unrecognized shard statuses"})
    complexes = read_tsv(analysis / "complexes" / "complex_inventory.tsv")
    failed_mapping = sum(r["status"] == "FAILED" for r in complexes)
    qc.append({"check": "ATOM_MAPPING", "status": "PASS" if complexes and failed_mapping == 0 else "FAIL", "observed": failed_mapping, "expected": 0, "message": "failed complex builds"})
    reps = read_tsv(analysis / "clustering" / "representatives.tsv")
    qc.append({"check": "GLOBAL_CLUSTERING", "status": "PASS" if reps else "FAIL", "observed": len(reps), "expected": ">0", "message": "representative rows"})
    plip = read_tsv(analysis / "interactions" / "plip_execution_status.tsv")
    failed_rep_plip = sum(r["kind"] == "representative" and r["status"] not in {"SUCCESS", "SKIPPED_COMPLETE"} for r in plip)
    qc.append({"check": "REPRESENTATIVE_PLIP", "status": "PASS" if plip and failed_rep_plip == 0 else "WARNING", "observed": failed_rep_plip, "expected": 0, "message": "failed/missing representative PLIP analyses; geometry and energy comparisons remain valid"})
    receptor_md5 = read_tsv(dock / "receptor_md5.tsv")
    md5s = defaultdict(set)
    for r in receptor_md5: md5s[(r["cas"], r["receptor"])].add(r["md5"])
    inconsistent = sum(len(v) > 1 for v in md5s.values())
    qc.append({"check": "RECEPTOR_MD5_CONSISTENCY", "status": "PASS" if inconsistent == 0 else "WARNING", "observed": inconsistent, "expected": 0, "message": "CAS/receptor groups with differing complete-receptor copies"})
    case_matches = read_tsv(dock / "case_insensitive_matches.tsv")
    qc.append({"check": "CASE_INSENSITIVE_AUDIT", "status": "PASS", "observed": len(case_matches), "expected": "AUDITED", "message": "automatic case-insensitive resolutions retained"})
    qc.extend(comparison_qc_rows(analysis))
    write_tsv(analysis / "quality_control.tsv", qc, ["check", "cas", "receptor", "pocket", "status", "observed", "expected", "message"])
    critical_fail = sum(r["status"] == "FAIL" for r in qc)
    update_stage_status(analysis, "validation", "PASS" if critical_fail == 0 else "FAIL", "{} QC failures".format(critical_fail), 0 if critical_fail == 0 else 1)
    return 0 if critical_fail == 0 else 1


if __name__ == "__main__": raise SystemExit(main())
