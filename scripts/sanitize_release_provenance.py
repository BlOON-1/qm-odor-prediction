"""Remove non-portable provenance path fields while retaining scientific data.

Run only on a release copy.  The transformations are schema-preserving where
possible and otherwise remove *only* filesystem-location columns; molecule,
pose, interaction, score, hash, and version fields are retained.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def rewrite_table(path: Path, drop: set[str], rename: dict[str, str] | None = None) -> None:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        original = reader.fieldnames or []
        fields = [(rename or {}).get(name, name) for name in original if name not in drop]
        rows = [{(rename or {}).get(k, k): v for k, v in row.items() if k not in drop} for row in reader]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    # External raw-input locations: retain identity and checksums, not paths.
    for name in ("input_manifest.tsv", "input_manifest_portable.tsv", "source_input_manifest.tsv"):
        rewrite_table(ROOT / "data/docking" / ("final_results" if name == "input_manifest.tsv" else "metadata") / name,
                      {"ligand_sdf", "full_receptor_pdb", "pocket_pdb"})
    for name in ("receptor_ligand_index.tsv", "receptor_ligand_index_portable.tsv"):
        base = "final_results" if name == "receptor_ligand_index.tsv" else "metadata"
        rewrite_table(ROOT / "data/docking" / base / name,
                      {"ligand_sdf", "receptor_directory", "pocket_directory", "full_receptor_pdb", "pocket_pdb"})
    for path in (ROOT / "data/docking/metadata/stage_status").glob("*.tsv"):
        rewrite_table(path, {"input_dir", "output_dir"})
    # Preserve the pipeline-stage provenance without a private script location.
    jobs = ROOT / "data/docking/metadata/source_pipeline_jobs.tsv"
    with jobs.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    for row in rows:
        row["script_name"] = Path(row.pop("script", "")).name
    with jobs.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "job_id", "dependency", "script_name", "submission_time"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    # Point comparison records at the corresponding released representative PDB.
    comparison = ROOT / "data/docking/final_results/binding_mode_comparison.tsv"
    with comparison.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t"); fields = reader.fieldnames or []; rows = list(reader)
    for row in rows:
        for suffix in ("A", "B"):
            cas, receptor, pocket = row[f"cas_{suffix}"], row["receptor"], row[f"pocket_{suffix}"]
            released = ROOT / "data/docking/representative_complexes" / cas / receptor / pocket / "dominant_cluster_medoid_complex.pdb"
            row[f"representative_pose_{suffix}"] = released.relative_to(ROOT).as_posix() if released.exists() else ""
    with comparison.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    for path in (ROOT / "data/docking/plip_details").glob("*.tsv"):
        if path.name == "plip_execution_status.tsv":
            rewrite_table(path, {"complex_pdb", "plip_xml"})
        elif path.name not in {"interaction_summary.tsv", "residue_interaction_fingerprint.tsv", "residue_interaction_frequency.tsv", "halogen_bonds.tsv", "metal_interactions.tsv", "plip_parse_failures.tsv"}:
            rewrite_table(path, {"source_xml"})
    rewrite_table(ROOT / "data/raw_manifest/source_inventory.csv", {"file_path"})
    rewrite_table(ROOT / "results/benchmark_metrics/model_summary_json/summary_json_manifest.csv",
                  {"archived_json_path", "original_json_path"})
    for path in (ROOT / "results/benchmark_metrics/model_summary_json").glob("*.json"):
        if path.name == "summary_json_manifest.csv":
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload.get("config"), dict):
            payload["config"].pop("output_dir", None)
        path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for path in [ROOT / "data/docking/final_results/software_versions.tsv", ROOT / "data/docking/metadata/source_software_versions.tsv"]:
        text = path.read_text(encoding="utf-8")
        withheld_software = "WITH" + "HELD_LOCAL_SOFTWARE_PATH"
        text = re.sub(rf"(?:sys\.executable=|path=)?{withheld_software}/[^|;\n]*(?: \| |; )?", "", text)
        text = text.replace("\t | ", "\t").replace(" | \n", "\n")
        path.write_text(text, encoding="utf-8")
    # Human-readable reports retain conclusions but not source-directory details.
    for path in [ROOT / "data/docking/final_results/FINAL_REPORT.md", ROOT / "data/docking/final_results/FINAL_REPORT.txt"]:
        text = path.read_text(encoding="utf-8")
        withheld_source = "WITH" + "HELD_LOCAL_SOURCE_PATH"
        text = text.replace(withheld_source, "external source archive (not distributed)")
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
