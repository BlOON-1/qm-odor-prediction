#!/usr/bin/env python3
"""Stage a portable, read-only release subset of the ESRS docking outputs.

This utility deliberately never writes below SOURCE_ROOT.  It copies only the
tabular result products, selected representative structures, and the scripts
that produced those products; large raw pose/XML directories are excluded.
"""
from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("DOCKING_SOURCE_ROOT", "NOT_CONFIGURED"))
SOURCE_SCRIPTS = Path(os.environ.get("DOCKING_SOURCE_SCRIPTS", "NOT_CONFIGURED"))
OUT = PROJECT / "data" / "docking"
REPORT = PROJECT / "reports" / "docking_data_audit_report.md"
DICTIONARY = PROJECT / "docs" / "docking_data_dictionary.md"

CORE = [
    "FINAL_REPORT.md", "FINAL_REPORT.txt", "PIPELINE_STATUS.tsv", "quality_control.tsv",
    "energy_summary.tsv", "cluster_summary.tsv", "cutoff_sensitivity.tsv",
    "interaction_summary.tsv", "residue_interaction_frequency.tsv",
    "binding_mode_comparison.tsv", "receptor_ligand_index.tsv",
    "receptor_ligand_groups.tsv", "comparison_eligibility.tsv", "comparison_skips.tsv",
    "comparison_discovery_status.tsv", "all_failures.tsv", "software_versions.tsv",
    "input_manifest.tsv", "output_manifest.tsv", "RESULT_FILES_README.txt",
]
PLIP = [
    "hydrogen_bonds.tsv", "hydrophobic_contacts.tsv", "salt_bridges.tsv", "pi_stacking.tsv",
    "pi_cation.tsv", "halogen_bonds.tsv", "metal_interactions.tsv",
    "residue_interaction_fingerprint.tsv", "plip_execution_status.tsv",
    "plip_parse_failures.tsv", "interaction_summary.tsv", "residue_interaction_frequency.tsv",
]
SCRIPTS = [
    "run_autodock4_batch_pipeline_64core.sh", "extract_dlg_poses.py", "cluster_poses.py",
    "build_complexes.py", "run_plip_batch.py", "parse_plip_results.py",
    "compare_binding_modes.py", "generate_final_report.py", "validate_pipeline_outputs.py",
    "pipeline_config.sh", "pipeline_common.py",
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path: Path, fieldnames, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(data)

def copy(source: Path, dest: Path, transfers: list):
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
    transfers.append((str(source), str(dest.relative_to(PROJECT)).replace("\\", "/"), source.stat().st_size, sha256(source), sha256(dest)))

def portable_value(value: str) -> str:
    # Existing source paths are provenance-only and their targets are not shipped.
    return "" if re.match(r"(?:[A-Za-z]:\\|/home/)", value or "") else value

def main():
    if not SOURCE_ROOT.exists(): raise SystemExit(f"Source unavailable: {SOURCE_ROOT}")
    if OUT.exists(): raise SystemExit(f"Refusing to overwrite existing release directory: {OUT}")
    transfers, missing = [], []
    final_dir, plip_dir, meta_dir = OUT / "final_results", OUT / "plip_details", OUT / "metadata"
    for name in CORE:
        src = SOURCE_ROOT / "final_results" / name
        if src.exists(): copy(src, final_dir / name, transfers)
        else: missing.append(f"final_results/{name}")
    for name in PLIP:
        src = SOURCE_ROOT / "interactions" / name
        if src.exists(): copy(src, plip_dir / name, transfers)
        else: missing.append(f"interactions/{name}")
    # Metadata is separately copied from the simulation root, retaining raw source bytes.
    for name in ("input_manifest.tsv", "software_versions.tsv", "pipeline_jobs.tsv"):
        src = SOURCE_ROOT / name
        if src.exists(): copy(src, meta_dir / f"source_{name}", transfers)
        else: missing.append(name)
    stage_dir = SOURCE_ROOT / "stage_status"
    if stage_dir.exists():
        for src in sorted(stage_dir.glob("*.tsv")): copy(src, meta_dir / "stage_status" / src.name, transfers)
    script_dir = PROJECT / "scripts" / "docking"
    for name in SCRIPTS:
        src = SOURCE_SCRIPTS / name
        if src.exists(): copy(src, script_dir / name, transfers)
        else: missing.append(f"pipeline_scripts/{name}")

    # Publish a clean index without unshipped, machine-specific absolute paths.
    source_index = SOURCE_ROOT / "final_results" / "receptor_ligand_index.tsv"
    index_rows = rows(source_index)
    public_index = []
    for r in index_rows:
        public_index.append({k: portable_value(v) for k, v in r.items()})
    write_tsv(meta_dir / "receptor_ligand_index_portable.tsv", list(index_rows[0]), public_index)
    source_manifest = SOURCE_ROOT / "input_manifest.tsv"
    if source_manifest.exists():
        manifest_rows = rows(source_manifest)
        write_tsv(meta_dir / "input_manifest_portable.tsv", list(manifest_rows[0]), [{k: portable_value(v) for k, v in r.items()} for r in manifest_rows])

    reps = rows(SOURCE_ROOT / "clustering" / "representatives.tsv")
    rep_index, structural_issues = [], []
    for r in reps:
        cas, receptor, pocket, role = r["cas"], r["receptor"], r["pocket"], r["role"]
        filename = f"{role}_complex.pdb"
        src = SOURCE_ROOT / "complexes" / cas / receptor / pocket / "representative" / filename
        dest = OUT / "representative_complexes" / cas / receptor / pocket / filename
        released = str(dest.relative_to(PROJECT)).replace("\\", "/")
        if src.exists():
            copy(src, dest, transfers)
            pdb = src.read_text(encoding="utf-8", errors="replace")
            receptor_atoms = sum(1 for x in pdb.splitlines() if x.startswith("ATOM"))
            ligand_atoms = sum(1 for x in pdb.splitlines() if x.startswith("HETATM"))
            if not receptor_atoms or not ligand_atoms: structural_issues.append(f"{released}: ATOM={receptor_atoms}, HETATM={ligand_atoms}")
            digest = sha256(dest)
        else:
            digest = ""; structural_issues.append(f"Missing representative: {src}")
        rep_index.append({"CAS": cas, "Receptor": receptor, "Pocket": pocket, "Pose_Type": role,
                          "Pose_ID": r["pose_id"], "Cluster_ID": "cluster_001" if role == "dominant_cluster_medoid" else "",
                          "Source_File": "WITHHELD_ABSOLUTE_SOURCE_PATH", "Released_File": released if src.exists() else "",
                          "SHA256": digest})
    write_tsv(meta_dir / "representative_pose_index.tsv", list(rep_index[0]), rep_index)

    # Validation products: report discoverability only, never invent/reconstruct them.
    validation_candidates = [p for p in SOURCE_ROOT.rglob("*") if p.is_file() and re.search(r"8F76|8UXY|9WG4", str(p), re.I)]
    # Tests and counts.
    qc, energy, cluster = rows(SOURCE_ROOT / "final_results" / "quality_control.tsv"), rows(SOURCE_ROOT / "final_results" / "energy_summary.tsv"), rows(SOURCE_ROOT / "final_results" / "cluster_summary.tsv")
    plip_status = rows(SOURCE_ROOT / "interactions" / "plip_execution_status.tsv")
    fingerprint = rows(SOURCE_ROOT / "interactions" / "residue_interaction_fingerprint.tsv")
    frequency = rows(SOURCE_ROOT / "interactions" / "residue_interaction_frequency.tsv")
    # The authoritative index uses explicit original-name column names, unlike
    # the downstream summary tables that shorten them to receptor/pocket.
    combos = {(r["cas"], r["receptor_name_original"], r["pocket_name_original"]) for r in index_rows}
    ligand_receptors = {(r["cas"], r["receptor_name_original"]) for r in index_rows}
    plip_success = [r for r in plip_status if r.get("status") == "SUCCESS"]
    interaction_a = SOURCE_ROOT / "final_results" / "interaction_summary.tsv"
    interaction_b = SOURCE_ROOT / "interactions" / "interaction_summary.tsv"
    frequency_a = SOURCE_ROOT / "final_results" / "residue_interaction_frequency.tsv"
    frequency_b = SOURCE_ROOT / "interactions" / "residue_interaction_frequency.tsv"
    duplicate_status = {"interaction_summary": sha256(interaction_a) == sha256(interaction_b), "residue_interaction_frequency": sha256(frequency_a) == sha256(frequency_b)}
    duplicate_combination_rows = len(index_rows) - len(combos)
    # Direct row duplication is meaningful only using complete records, not a selected key.
    dup_counts = {"quality_control": len(qc) - len({tuple(sorted(x.items())) for x in qc}), "cluster_summary": len(cluster) - len({tuple(sorted(x.items())) for x in cluster})}
    # Figure 4 ligand-to-receptor mapping comes from numeric CAS/receptor association; names are not in this source index.
    expected = {"OR1D2": {"106-25-2", "141-12-8", "18127-01-0"}, "OR1A1": {"66-25-1", "124-19-6", "63767-86-2"}, "OR2W1": {"71-41-0", "111-29-5", "20602-31-7"}}
    observed = defaultdict(set)
    for cas, receptor in ligand_receptors: observed[receptor].add(cas)
    mapping_checks = {r: {"expected": sorted(v), "observed": sorted(observed[r]), "match": v == observed[r]} for r, v in expected.items()}
    # Manifest of every copied immutable source object.
    write_tsv(meta_dir / "file_transfer_manifest.tsv", ["source_path_withheld", "released_file", "bytes", "source_sha256", "released_sha256"],
              [{"source_path_withheld": "WITHHELD_ABSOLUTE_SOURCE_PATH", "released_file": b, "bytes": c, "source_sha256": d, "released_sha256": e} for a,b,c,d,e in transfers])
    # File-level TSV inventories (actual headers/rows) support independent inspection.
    inventory = []
    for p in sorted(OUT.rglob("*.tsv")):
        try:
            rr = rows(p); hdr = list(rr[0]) if rr else (p.read_text(encoding="utf-8-sig").splitlines()[0].split("\t") if p.stat().st_size else [])
            inventory.append({"file": str(p.relative_to(PROJECT)).replace("\\", "/"), "bytes": p.stat().st_size, "data_rows": len(rr), "header": "|".join(hdr), "sha256": sha256(p)})
        except Exception as exc: inventory.append({"file": str(p.relative_to(PROJECT)), "bytes": p.stat().st_size, "data_rows": "READ_ERROR", "header": str(exc), "sha256": sha256(p)})
    write_tsv(meta_dir / "released_tsv_inventory.tsv", list(inventory[0]), inventory)

    README = f"""# Docking and PLIP release subset

## Overview

This directory is a portable, review-oriented subset staged from an immutable local AutoDock4/PLIP run. It contains final tabular outputs, pose-level PLIP details, {len(rep_index)} representative receptor--ligand complexes, provenance hashes, and the analysis scripts used by the source run. It intentionally excludes raw PLIP XML, all-pose PDBs, DLG/AD4 runs, logs, and other large intermediates.

## Systems and scope

The source index contains {len({r['cas'] for r in index_rows})} CAS identifiers, {len({r['receptor_name_original'] for r in index_rows})} receptors, {len(combos)} receptor--ligand--pocket combinations, and {len(ligand_receptors)} receptor--ligand pairs. Each indexed pair is evaluated in four pockets. The specific Figure-4 name/CAS mapping needs human confirmation from a ligand-name source; this package preserves the authoritative CAS/receptor mappings.

## Contents

- `final_results/`: report, QC, energy, clustering, comparison, and summary tables copied byte-for-byte from the source final-results directory.
- `plip_details/`: interaction event tables, pose fingerprints, PLIP execution status, parse failures, and duplicate summary tables.
- `representative_complexes/`: global-lowest-energy and dominant-cluster-medoid PDBs. Coordinates were copied without modification.
- `metadata/`: raw-source metadata copies, portable indexes, per-file SHA-256 transfer manifest, TSV inventory, and stage status.
- `../../scripts/docking/`: archived pipeline scripts, copied byte-for-byte.

## Portable paths and provenance

`metadata/*_portable.tsv` clears non-portable absolute source-path fields; it is a derived release index, not byte-identical to its source. Original source file hashes and released hashes are in `metadata/file_transfer_manifest.tsv`. The core `final_results` copies retain their source bytes, including historical paths; do not treat those locations as usable paths.

## Representative structures

Selection keys are in `metadata/representative_pose_index.tsv`. `global_lowest_energy` selects the globally lowest extracted pose; `dominant_cluster_medoid` selects the medoid from the dominant cluster. A docking score or PLIP contact is computational evidence, not experimental receptor binding or activation evidence.

## Reproduction and boundaries

Run `python scripts/audit_and_stage_docking.py` only to stage a fresh package from the original local source. It does not run AutoDock4 or PLIP. For Figure 4F review, link `plip_details/residue_interaction_fingerprint.tsv` and `plip_details/residue_interaction_frequency.tsv` by CAS/receptor/pocket/pose_id, then apply the inclusion criteria in `final_results/binding_mode_comparison.tsv` and `comparison_eligibility.tsv`. See the audit report for observed limitations.
"""
    (OUT / "README.md").write_text(README, encoding="utf-8")
    DICTIONARY.write_text(f"""# Docking data dictionary

## Common keys

`cas`, `receptor`, and `pocket` identify an indexed docking system; `pose_id` identifies an extracted pose within that system. `kind` and `role` distinguish all-pose versus representative/selection contexts. Missing strings are source-provided empty values unless stated otherwise.

## Principal tables

- `final_results/receptor_ligand_index.tsv`: input identity and pocket/grid metadata. `gridcenter`, `npts`, and `spacing` define the AutoDock grid; source path fields are historical and non-portable.
- `final_results/energy_summary.tsv`: one row per CAS/receptor/pocket; `pose_count`, minimum, median, mean and SD docking energies (kcal/mol).
- `final_results/cluster_summary.tsv`: one row per cluster; `cluster_fraction` is cluster pose count divided by indexed pose count; medoid and lowest-energy pose IDs are explicit.
- `final_results/quality_control.tsv`: individual QC assertions with observed and expected values.
- `final_results/binding_mode_comparison.tsv`, `comparison_eligibility.tsv`, and `comparison_skips.tsv`: cross-ligand comparison outputs and inclusion/exclusion accounting.
- `plip_details/*_contacts.tsv` and `*_bonds.tsv`: one PLIP interaction event per row.
- `plip_details/residue_interaction_fingerprint.tsv`: one binary residue/type feature per pose (`present` = 1).
- `plip_details/residue_interaction_frequency.tsv`: a residue/type aggregation. `pose_occurrence_count` is divided by `frequency_denominator`; requested and successful PLIP pose counts are retained so the denominator is auditable.
- `plip_details/plip_execution_status.tsv`: one PLIP invocation/status record per complex; `SUCCESS` means XML generation, not experimental validation.
- `metadata/representative_pose_index.tsv`: publication-facing mapping from selection type and pose ID to copied PDB, including release SHA-256.
- `metadata/file_transfer_manifest.tsv`: byte-level source/release identity check for copied files. Source locations are deliberately withheld in the public release.

## File relationships

`residue_interaction_fingerprint` and interaction-event tables join to `plip_execution_status` using CAS/receptor/pocket/pose_id. Representative PDBs join through `representative_pose_index`. For cross-ligand figures, use comparison eligibility and the same pocket/selection scope; do not combine all-pose frequencies with medoid-only data without explicitly changing the denominator.
""", encoding="utf-8")
    report = f"""# Docking/PLIP data audit report

## Source and staging

Source examined: local ESRS docking simulation and pipeline-script directories (not altered). Release staged at `data/docking/`. {len(transfers)} files were copied and checked with SHA-256; all copied source/release hash pairs match: **{all(d == e for a,b,c,d,e in transfers)}**. Missing requested files: {', '.join(missing) if missing else 'none'}.

## Observed scale

- Indexed CAS values: **{len({r['cas'] for r in index_rows})}**; receptors: **{len({r['receptor_name_original'] for r in index_rows})}**; pockets: **{len({r['pocket_name_original'] for r in index_rows})}**.
- Receptor--ligand pairs: **{len(ligand_receptors)}**; receptor--ligand--pocket combinations: **{len(combos)}**. The index has {duplicate_combination_rows} duplicate complete-combination rows.
- Energy-summary rows: **{len(energy)}**; QC rows: **{len(qc)}**; cluster-summary rows: **{len(cluster)}**.
- PLIP execution records: **{len(plip_status)}**; `SUCCESS`: **{len(plip_success)}**; non-success: **{len(plip_status)-len(plip_success)}**. Fingerprint rows: **{len(fingerprint)}**; residue-frequency rows: **{len(frequency)}**.
- Representative index rows/PDBs expected: **{len(rep_index)}**; missing or invalid receptor/ligand PDB checks: **{len(structural_issues)}**.

## Consistency checks

- Same-byte duplicate: final-results vs interactions `interaction_summary.tsv`: **{duplicate_status['interaction_summary']}**; `residue_interaction_frequency.tsv`: **{duplicate_status['residue_interaction_frequency']}**. The release keeps both source copies in their logical folders and identifies their duplicate status here; downstream publication may retain one canonical copy only after review.
- Exact duplicate rows: quality-control={dup_counts['quality_control']}; cluster-summary={dup_counts['cluster_summary']}.
- Expected Figure-4 CAS/receptor sets inferred from the supplied task: `{mapping_checks}`. These are CAS-only checks; molecule-name/CAS identity and experimental-evidence status are **not present in the inspected docking index** and require confirmation against an authoritative ligand/evidence table.
- Validation-file discovery for PDB 8F76/8UXY/9WG4 within `docking_simulation`: **{len(validation_candidates)}** candidate files. No validation files were copied automatically; a candidate must contain usable alignment/RMSD/pose-recovery data before inclusion.

## Interpretation and release boundaries

The results support audit of docking score summaries, clustering, PLIP event/fingerprint/frequency records, and binding-mode comparison inputs. They do **not** alone establish experimental receptor recognition/activation or successful experimental-pose recovery. In particular, this staging run does not label the pipeline as a successful structural validation. Raw XML, all-pose PDBs, AutoDock run files, and logs remain at the source location; archive those separately (for example, checksummed split archives in Zenodo) if full rerun provenance is required.

## Follow-up needed

1. Confirm molecule names, CAS assignments, literature receptor evidence, and Figure-4 panel pose selection from the manuscript/source-data table.
2. Inspect the {len(validation_candidates)} validation candidates for 8F76/8UXY/9WG4 before creating a public `validation/` subset.
3. Review licensing/redistribution terms for input receptor/ligand structures and software before public release.
4. Treat original absolute paths as confidential/non-portable; use `metadata/*_portable.tsv` for public processing.
"""
    REPORT.write_text(report, encoding="utf-8")
    print(f"Staged {len(transfers)} files; {len(rep_index)} representative PDBs; report={REPORT}")

if __name__ == "__main__": main()
