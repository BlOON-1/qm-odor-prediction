#!/usr/bin/env python3
"""Generate a traceable final report even when an upstream stage failed."""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np

from pipeline_common import analysis_run_dir, dock_run_dir, read_tsv, update_stage_status, write_tsv


def copy_or_empty(source, target, fields):
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_file() and source.stat().st_size:
        shutil.copy2(str(source), str(target))
    else:
        write_tsv(target, [], fields)


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--analysis-run", type=Path, default=None); parser.add_argument("--dock-run", type=Path, default=None); args = parser.parse_args()
    analysis, dock = args.analysis_run or analysis_run_dir(), args.dock_run or dock_run_dir()
    final = analysis / "final_results"; final.mkdir(parents=True, exist_ok=True)

    poses = [r for r in read_tsv(analysis / "docking_post_analysis" / "pose_inventory.tsv") if r.get("extraction_status") == "SUCCESS"]
    energy_groups = {}
    for row in poses:
        key = (row["cas"], row["receptor"], row["pocket"]); energy_groups.setdefault(key, []).append(float(row["binding_energy_kcal_mol"]))
    energy_summary = []
    for key, values in sorted(energy_groups.items()):
        energy_summary.append({"cas": key[0], "receptor": key[1], "pocket": key[2], "pose_count": len(values),
                               "minimum_kcal_mol": "{:.3f}".format(min(values)), "median_kcal_mol": "{:.3f}".format(float(np.median(values))),
                               "mean_kcal_mol": "{:.3f}".format(float(np.mean(values))), "sd_kcal_mol": "{:.3f}".format(float(np.std(values, ddof=1)) if len(values) > 1 else 0.0)})
    write_tsv(final / "energy_summary.tsv", energy_summary, ["cas", "receptor", "pocket", "pose_count", "minimum_kcal_mol", "median_kcal_mol", "mean_kcal_mol", "sd_kcal_mol"])

    copies = [
        (analysis / "quality_control.tsv", "quality_control.tsv", ["check", "status", "message"]),
        (analysis / "clustering" / "cluster_summary.tsv", "cluster_summary.tsv", ["cas", "receptor", "pocket", "cluster_id"]),
        (analysis / "clustering" / "cutoff_sensitivity.tsv", "cutoff_sensitivity.tsv", ["cas", "receptor", "pocket", "cutoff_A"]),
        (analysis / "interactions" / "interaction_summary.tsv", "interaction_summary.tsv", ["cas", "receptor", "pocket", "interaction_type"]),
        (analysis / "interactions" / "residue_interaction_frequency.tsv", "residue_interaction_frequency.tsv", ["cas", "receptor", "pocket", "frequency"]),
        (analysis / "comparisons" / "receptor_ligand_index.tsv", "receptor_ligand_index.tsv", ["cas", "receptor_name_original", "canonical_receptor_key", "pocket_index", "input_status"]),
        (analysis / "comparisons" / "receptor_ligand_groups.tsv", "receptor_ligand_groups.tsv", ["canonical_receptor_key", "display_receptor_name", "unique_cas_count", "status"]),
        (analysis / "comparisons" / "comparison_eligibility.tsv", "comparison_eligibility.tsv", ["comparison_id", "canonical_receptor_key", "cas_A", "cas_B", "pocket_index", "eligible", "status"]),
        (analysis / "comparisons" / "comparison_discovery_status.tsv", "comparison_discovery_status.tsv", ["status", "receptor_group_count", "eligible_pocket_comparison_count", "message"]),
        (analysis / "comparisons" / "binding_mode_comparison.tsv", "binding_mode_comparison.tsv", ["comparison_id", "cas_A", "cas_B", "receptor", "pocket_index", "comparison_status"]),
        (analysis / "comparisons" / "comparison_skips.tsv", "comparison_skips.tsv", ["canonical_receptor_key", "cas_A", "cas_B", "pocket_index", "status", "reason"]),
        (analysis / "software_versions.tsv", "software_versions.tsv", ["component", "status", "version_or_message"]),
        (analysis / "input_manifest.tsv", "input_manifest.tsv", ["cas", "receptor", "pocket", "ligand_sdf", "full_receptor_pdb", "pocket_pdb"]),
    ]
    for source, name, fields in copies: copy_or_empty(source, final / name, fields)

    failures = []
    failure_sources = [dock / "failures.tsv", analysis / "docking_post_analysis" / "pose_extraction_failures.tsv",
                       analysis / "clustering" / "clustering_failures.tsv", analysis / "interactions" / "plip_parse_failures.tsv"]
    for source in failure_sources:
        for row in read_tsv(source): failures.append({"source": str(source), "record": " | ".join("{}={}".format(k, v) for k, v in row.items())})
    for row in read_tsv(analysis / "complexes" / "complex_inventory.tsv"):
        if row.get("status") == "FAILED": failures.append({"source": "complex_inventory.tsv", "record": "{}|{}|{}|{}".format(row["cas"], row["receptor"], row["pocket"], row["message"])})
    for row in read_tsv(analysis / "interactions" / "plip_execution_status.tsv"):
        if row.get("status") == "FAILED": failures.append({"source": "plip_execution_status.tsv", "record": "{}|{}|{}|{}".format(row["cas"], row["receptor"], row["pocket"], row["message"])})
    write_tsv(final / "all_failures.tsv", failures, ["source", "record"])

    stages = []
    for path in sorted((analysis / "stage_status").glob("*.tsv")):
        stages.extend(read_tsv(path))
    # Ignore a prior report-generation status when regenerating this report.
    stages = [row for row in stages if row.get("stage") != "report"]
    qc = read_tsv(analysis / "quality_control.tsv")
    required_stages = {"docking", "post_analysis", "extraction", "clustering", "complexes", "plip", "receptor_index", "comparison", "validation"}
    present_stages = {row.get("stage") for row in stages}
    originally_missing = sorted(required_stages - present_stages)
    for missing in sorted(required_stages - present_stages):
        stages.append({"stage": missing, "status": "MISSING", "exit_code": 1,
                       "message": "required stage status is missing; stage failed or never ran"})
    if not qc:
        qc = [{"check": "QUALITY_CONTROL_TABLE", "status": "FAIL", "message": "quality_control.tsv is missing or empty"}]
        write_tsv(final / "quality_control.tsv", qc,
                  ["check", "cas", "receptor", "pocket", "status", "observed", "expected", "message"])
    stage_by_name = {row.get("stage"): row.get("status", "MISSING") for row in stages}
    if originally_missing or stage_by_name.get("docking") == "MISSING" or stage_by_name.get("post_analysis") == "MISSING":
        pipeline_outcome = "STATUS_MISSING"
    elif stage_by_name.get("docking") != "SUCCESS":
        pipeline_outcome = "DOCKING_FAILED"
    elif stage_by_name.get("post_analysis") != "SUCCESS":
        pipeline_outcome = "POST_ANALYSIS_FAILED"
    elif any(stage_by_name.get(name) not in {"SUCCESS", "SUCCESS_NO_ELIGIBLE_PAIRS", "PASS", "WARNING"} for name in required_stages):
        pipeline_outcome = "POST_ANALYSIS_FAILED"
    else:
        pipeline_outcome = "ALL_SUCCESS"
    failure_statuses = {"FAIL", "FAILED", "MISSING", "SKIPPED", "RUNNING"}
    overall = "FAIL" if pipeline_outcome != "ALL_SUCCESS" or any(r.get("status") in failure_statuses for r in qc + stages) else ("WARNING" if failures or any(r.get("status") == "WARNING" for r in qc + stages) else "PASS")
    pipeline_status = [{"stage": r.get("stage", "UNKNOWN"), "status": r.get("status", "UNKNOWN"), "exit_code": r.get("exit_code", ""), "message": r.get("message", "")} for r in stages]
    pipeline_status.append({"stage": "report", "status": "PASS", "exit_code": 0,
                            "message": "final report artifacts were generated; pipeline outcome is recorded separately"})
    pipeline_status.append({"stage": "PIPELINE_OUTCOME", "status": pipeline_outcome,
                            "exit_code": 0 if pipeline_outcome == "ALL_SUCCESS" else 1,
                            "message": "distinguishes docking failure, post-analysis failure, complete success, and missing status"})
    pipeline_status.append({"stage": "OVERALL", "status": overall, "exit_code": 0 if overall in {"PASS", "WARNING"} else 1, "message": "derived from stage status, QC, and failure tables"})
    write_tsv(final / "PIPELINE_STATUS.tsv", pipeline_status, ["stage", "status", "exit_code", "message"])

    manifests = read_tsv(dock / "job_manifest.tsv")
    cluster_rows = read_tsv(analysis / "clustering" / "cluster_summary.tsv")
    sensitivity = read_tsv(analysis / "clustering" / "cutoff_sensitivity.tsv")
    plip_rows = read_tsv(analysis / "interactions" / "plip_execution_status.tsv")
    plip_success = sum(r.get("status") in {"SUCCESS", "SKIPPED_COMPLETE"} for r in plip_rows)
    dominant = [r for r in cluster_rows if r.get("cluster_id") == "cluster_001"]
    frequencies = sorted(read_tsv(analysis / "interactions" / "residue_interaction_frequency.tsv"), key=lambda r: float(r.get("frequency") or 0), reverse=True)[:20]
    comparisons = read_tsv(analysis / "comparisons" / "binding_mode_comparison.tsv")
    receptor_groups = read_tsv(analysis / "comparisons" / "receptor_ligand_groups.tsv")
    comparison_eligibility = read_tsv(analysis / "comparisons" / "comparison_eligibility.tsv")
    comparison_skips = read_tsv(analysis / "comparisons" / "comparison_skips.tsv")
    single_receptor_groups = [row for row in receptor_groups if row.get("status") == "SKIPPED_SINGLE_LIGAND"]
    multi_receptor_groups = [row for row in receptor_groups if int(row.get("unique_cas_count") or 0) >= 2]
    eligible_receptor_groups = [row for row in receptor_groups if row.get("comparison_eligible") == "1"]
    theoretical_pairs = sum(int(row.get("theoretical_ligand_pair_count") or 0) for row in receptor_groups)
    common_pocket_rows = [row for row in comparison_eligibility if row.get("pocket_A") and row.get("pocket_B")]
    complete_comparisons = [row for row in comparisons if row.get("comparison_status") == "SUCCESS"]
    partial_comparisons = [row for row in comparisons if row.get("comparison_status", "").startswith("PARTIAL")]
    skip_reasons = Counter(row.get("status") or "UNKNOWN" for row in comparison_skips)
    qm_file = analysis / "comparisons" / "qm_analysis_status.txt"
    qm_status = qm_file.read_text(errors="replace").strip() if qm_file.is_file() else "QM_ANALYSIS_STATUS=NOT_PROVIDED"
    params = {r["parameter"]: r["value"] for r in read_tsv(dock / "run_parameters.tsv")}

    lines = ["# FINAL REPORT — AutoDock4 High-Precision Batch Pipeline", "", "## 1. Run information", "",
             "- Pipeline ID: `{}`".format(analysis.name), "- Generated: {}".format(datetime.now().astimezone().isoformat()),
             "- Pipeline execution outcome: **{}**".format(pipeline_outcome),
             "- Overall quality status: **{}**".format(overall), "- Dock result directory: `{}`".format(dock), "",
             "## 2. Software and input scope", "", "- Configured CAS–receptor–pocket combinations: **{}**".format(len(manifests)),
             "- CAS count: **{}**; receptor count: **{}**; pocket count: **{}**".format(len({r.get('cas') for r in manifests}), len({r.get('receptor') for r in manifests}), len({(r.get('receptor'), r.get('pocket')) for r in manifests})),
             "- Exact versions are recorded in `software_versions.tsv`.", "", "## 3. Docking parameters and CPU allocation", "",
             "- AutoDock workers: **{}**; post-analysis CPUs: **32**; PLIP worker cap: **16**.".format(params.get("max_parallel", "64")),
             "- LGA runs/pocket: **{}**; evaluations/run: **{}**; generations/run: **{}**; population: **{}**; RMSTOL: **{} Å**.".format(params.get("ga_run_total_per_pocket", "100"), params.get("ga_num_evals_per_run", "10000000"), params.get("ga_num_generations_per_run", "50000"), params.get("ga_population_size", "300"), params.get("rmstol_A", "2.0")),
             "", "## 4. Completion and energy statistics", "", "- Validated extracted poses: **{}** across **{}** groups.".format(len(poses), len(energy_groups)),
             "- Per-pocket 100-search completion and energy distributions are in `energy_summary.tsv` and `quality_control.tsv`.", "",
             "## 5. Global pose clustering", "", "- Primary cutoff: **2.0 Å**, fixed receptor coordinate frame, average linkage.",
             "- Dominant clusters summarized: **{}**. Sensitivity rows (1.5/2.0/2.5 Å): **{}**.".format(len(dominant), len(sensitivity)),
             "- Lowest-energy and dominant-cluster medoid identity is traceable via `representatives.tsv` and `pose_cluster_assignments.tsv`.", "",
             "## 6. PLIP interactions", "", "- Successful/reused PLIP analyses: **{}/{}**.".format(plip_success, len(plip_rows)),
             "- Frequencies use the actual successful all-pose PLIP count as denominator, never an assumed 100.", "- High-frequency residue rows reported: **{}** (top rows shown below).".format(len(frequencies)), ""]
    for row in frequencies[:10]:
        lines.append("- {} / {} / {}: {} {}{}{} — frequency {} (denominator {}).".format(row.get("cas"), row.get("receptor"), row.get("pocket"), row.get("interaction_type"), row.get("residue_chain"), row.get("residue_name"), row.get("residue_number"), row.get("frequency"), row.get("frequency_denominator")))
    lines += ["", "PLIP salt bridges, AutoDock4 electrostatic score components, and simple charged-residue proximity are distinct concepts; this workflow does not relabel proximity as an electrostatic interaction.", "",
              "## 7. Cross-ligand comparison discovery", "",
              "- Receptor names discovered: **{}**.".format(len(receptor_groups)),
              "- Single-CAS receptor groups: **{}**.".format(len(single_receptor_groups)),
              "- Receptor groups with two or more CAS directories: **{}**.".format(len(multi_receptor_groups)),
              "- Comparison-eligible receptor groups: **{}**.".format(len(eligible_receptor_groups)),
              "- Theoretical unique CAS pairs: **{}**.".format(theoretical_pairs),
              "- Actual common-pocket eligibility rows: **{}**.".format(len(common_pocket_rows)),
              "- Skipped records: **{}**.".format(len(comparison_skips)),
              "- Complete comparisons: **{}**; partial-data comparisons: **{}**.".format(len(complete_comparisons), len(partial_comparisons)),
              "- Skip reasons: {}.".format("; ".join("{}={}".format(key, value) for key, value in sorted(skip_reasons.items())) or "none"), ""]
    if not multi_receptor_groups:
        lines += ["No receptor name was shared by two or more unique CAS directories. Cross-ligand comparison was not applicable for this input set.", ""]
    else:
        lines += ["| Receptor | CAS count | CAS list | Ligand pairs | Common pockets | Successful | Skipped |",
                  "|---|---:|---|---:|---:|---:|---:|"]
        for group in multi_receptor_groups:
            canonical = group.get("canonical_receptor_key")
            related_eligibility = [row for row in comparison_eligibility if row.get("canonical_receptor_key") == canonical]
            related_results = [row for row in comparisons if row.get("canonical_receptor_key") == canonical]
            related_skips = [row for row in comparison_skips if row.get("canonical_receptor_key") == canonical]
            lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                group.get("display_receptor_name"), group.get("unique_cas_count"), group.get("cas_list"),
                group.get("theoretical_ligand_pair_count"),
                sum(bool(row.get("pocket_A") and row.get("pocket_B")) for row in related_eligibility),
                sum(row.get("comparison_status") == "SUCCESS" for row in related_results), len(related_skips)))
        lines.append("")
    lines += ["",
              "## 8. Cross-ligand binding-mode comparison", "", "- Completed/partial comparison rows: **{}**.".format(len(comparisons)),
              "- Metrics use centroid separation, principal-axis angle, residue overlap, interaction fingerprints, energy distributions, and dominant-cluster fractions—not ordinary cross-ligand all-atom RMSD.", "",
              "- Spatial-region labels are descriptive workflow classifications using parameters centroid <= {} A, residue Jaccard >= {}, and fingerprint similarity >= {}; they are not universal biological laws.".format(os.environ.get("COMPARISON_CENTROID_THRESHOLD_A", "5.0"), os.environ.get("COMPARISON_RESIDUE_JACCARD_THRESHOLD", "0.3"), os.environ.get("COMPARISON_FINGERPRINT_THRESHOLD", "0.3")), "",
              "## 9. Optional QM boundary", "", "- `{}`".format(qm_status), "- Missing QM inputs do not fail docking. No QM descriptors are invented from molecular structures.", "",
              "## 10. Failures, warnings, and QC", "", "- Recorded failure/warning records: **{}**.".format(len(failures)), "- See `PIPELINE_STATUS.tsv`, `quality_control.tsv`, and `all_failures.tsv`; an upstream critical failure prevents an all-success conclusion.", "",
              "## 11. Scientific interpretation boundaries", "",
              "1. AutoDock4 scores are model predictions, not experimental binding free energies.",
              "2. Molecular docking cannot prove receptor activation.",
              "3. A single lowest-energy pose does not establish the true binding mode.",
              "4. Dominant clusters and frequent interactions assess search reproducibility.",
              "5. Different ligands are not compared using ordinary all-atom RMSD.",
              "6. Docking supplies receptor-level structural-support evidence only.",
              "7. This workflow alone cannot establish that QM features cause odor differences."]
    (final / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    plain = [line.replace("**", "").replace("`", "").lstrip("# ") for line in lines]
    (final / "FINAL_REPORT.txt").write_text("\n".join(plain) + "\n", encoding="utf-8")

    readme = """FINAL RESULT FILES\n\nFINAL_REPORT.md / .txt: human-readable report.\nPIPELINE_STATUS.tsv and quality_control.tsv: completion truth; inspect first.\nenergy_summary.tsv: validated pose energy distributions.\ncluster_summary.tsv and cutoff_sensitivity.tsv: global clustering evidence.\ninteraction_summary.tsv and residue_interaction_frequency.tsv: PLIP evidence.\nbinding_mode_comparison.tsv: eligible cross-CAS comparisons within the same canonical receptor and pocket index.\nreceptor_ligand_index.tsv / receptor_ligand_groups.tsv: case-insensitive receptor discovery and reverse index.\ncomparison_eligibility.tsv / comparison_skips.tsv: structural/grid eligibility and traceable skips.\nall_failures.tsv: retained failures and warnings.\nsoftware_versions.tsv / input_manifest.tsv / output_manifest.tsv: provenance.\n"""
    (final / "RESULT_FILES_README.txt").write_text(readme, encoding="utf-8")
    output_rows = []
    for path in sorted(p for p in final.iterdir() if p.is_file() and p.name != "output_manifest.tsv"):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        output_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": digest})
    write_tsv(final / "output_manifest.tsv", output_rows, ["file", "bytes", "sha256"])
    update_stage_status(analysis, "report", "PASS", "Final report generated; pipeline overall status is {}".format(overall))
    return 0


if __name__ == "__main__": raise SystemExit(main())
