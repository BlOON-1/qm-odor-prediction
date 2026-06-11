"""Entry point for receptor annotation preparation and optional distance comparison."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qmodor.descriptors.common import backup_existing
from qmodor.receptor.compute_receptor_descriptor_distances import compute_receptor_descriptor_distances
from qmodor.receptor.prepare_m2or_receptor_annotations import prepare_m2or_receptor_annotations


def _upsert_yaml_files(project_root: Path) -> None:
    import yaml

    dataset_path = project_root / "configs" / "dataset.yaml"
    dataset = yaml.safe_load(dataset_path.read_text(encoding="utf-8")) if dataset_path.exists() else {}
    dataset.setdefault("files", {})
    dataset["files"].update(
        {
            "m2or_pair_or_raw": "data/receptor_annotations/m2or_pair_or_raw.csv",
            "m2or_matched_molecules": "data/receptor_annotations/m2or_matched_molecules.csv",
            "receptor_response_long": "data/receptor_annotations/receptor_response_long.csv",
            "receptor_response_vectors": "data/receptor_annotations/receptor_response_vectors.csv",
            "receptor_pairwise_concordance": "data/receptor_annotations/receptor_pairwise_concordance.csv",
            "receptor_concordant_divergent_pairs": "data/receptor_annotations/receptor_concordant_divergent_pairs.csv",
            "receptor_matching_exclusion_log": "data/receptor_annotations/receptor_matching_exclusion_log.csv",
            "structure_disjoint_train": "data/splits/structure_disjoint_splits/structure_disjoint_train.csv",
            "structure_disjoint_test": "data/splits/structure_disjoint_splits/structure_disjoint_test.csv",
            "structure_disjoint_split_audit": "data/splits/structure_disjoint_splits/structure_disjoint_split_audit.csv",
        }
    )
    backup_existing(dataset_path)
    dataset_path.write_text(yaml.safe_dump(dataset, sort_keys=False, allow_unicode=True), encoding="utf-8")

    figure_path = project_root / "configs" / "figure_settings.yaml"
    figure = yaml.safe_load(figure_path.read_text(encoding="utf-8")) if figure_path.exists() else {}
    figure.setdefault("figures", {})
    figure["figures"].setdefault("fig4", {})
    figure["figures"]["fig4"]["source_data"] = "data/source_data/fig4_source_data.xlsx"
    backup_existing(figure_path)
    figure_path.write_text(yaml.safe_dump(figure, sort_keys=False, allow_unicode=True), encoding="utf-8")


def _upsert_markdown(path: Path, heading: str, body: str) -> None:
    import re

    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    section = f"\n## {heading}\n\n{body.strip()}\n"
    pattern = re.compile(rf"\n## {re.escape(heading)}\n.*?(?=\n## |\Z)", re.S)
    if existing:
        new_text = pattern.sub(section, existing) if pattern.search(existing) else existing.rstrip() + section
    else:
        new_text = f"# {path.stem.replace('_', ' ').title()}\n" + section
    backup_existing(path)
    path.write_text(new_text.rstrip() + "\n", encoding="utf-8")


def _update_docs(project_root: Path) -> None:
    data_dict_path = project_root / "docs" / "data_dictionary.md"
    body = "\n".join(
        [
            "| File | Key fields or structure | Purpose |",
            "| --- | --- | --- |",
            "| `data/receptor_annotations/m2or_pair_or_raw.csv` | raw M2OR pairwise comparison rows | Local clean copy of the uploaded receptor-comparison source table. |",
            "| `data/receptor_annotations/m2or_matched_molecules.csv` | `cas`, `molecule_role_source`, count summaries | Molecule-level summary of matched receptor annotation coverage under the default species filter. |",
            "| `data/receptor_annotations/receptor_response_long.csv` | `pair_id`, `or_gene`, tested flags, response flags | Harmonized pair-receptor annotation table used for final comparison logic. |",
            "| `data/receptor_annotations/receptor_response_vectors.csv` | one row per CAS with OR-gene columns | Molecule-by-receptor response matrix for the default species filter. |",
            "| `data/receptor_annotations/receptor_pairwise_concordance.csv` | pair-level jointly tested receptor summary columns | Pairwise receptor concordance table used for orthogonal consistency grouping. |",
            "| `data/receptor_annotations/receptor_concordant_divergent_pairs.csv` | concordant/divergent pairs only | Final grouped pair table retained for group comparison. |",
            "| `data/receptor_annotations/receptor_matching_exclusion_log.csv` | exclusion reason log | Explicit log of rows excluded during duplicate handling, human-only filtering, or missing-data filtering. |",
            "| `data/splits/structure_disjoint_splits/structure_disjoint_train.csv` | split table with `cas` | Canonical structure-disjoint training split organized from uploaded fold files. |",
            "| `data/splits/structure_disjoint_splits/structure_disjoint_test.csv` | split table with `cas` | Canonical structure-disjoint test split organized from uploaded fold files. |",
            "| `data/splits/structure_disjoint_splits/structure_disjoint_split_audit.csv` | audit metrics | Leakage and role-convention audit for the uploaded structure-disjoint split files. |",
            "| `data/source_data/fig4_source_data.xlsx` | multiple sheets | Partial Fig. 4 source-data workbook using receptor annotation outputs and optional descriptor-distance comparisons. |",
        ]
    )
    _upsert_markdown(data_dict_path, "Receptor And Structure-Disjoint Files", body)


def _build_fig4_workbook(project_root: Path, receptor_result: dict, distance_result: dict) -> Path:
    import json

    source_data_dir = project_root / "data" / "source_data"
    source_data_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = source_data_dir / "fig4_source_data.xlsx"

    pairwise = pd.read_csv(project_root / "data" / "receptor_annotations" / "receptor_pairwise_concordance.csv")
    grouped = pd.read_csv(project_root / "data" / "receptor_annotations" / "receptor_concordant_divergent_pairs.csv")
    vectors = pd.read_csv(project_root / "data" / "receptor_annotations" / "receptor_response_vectors.csv")
    stats = pd.read_csv(project_root / "results" / "receptor_consistency" / "receptor_consistency_statistics.csv")
    if distance_result["generated"]:
        distance_sheet = pd.read_csv(project_root / "results" / "receptor_consistency" / "qm_rdkit_distance_by_receptor_group.csv")
    else:
        distance_sheet = pd.DataFrame(
            [
                {
                    "required_file": item,
                    "exists": False,
                    "status": "blocked",
                    "note": "Descriptor-distance comparison not generated.",
                }
                for item in distance_result["missing"]
            ]
        )
    audit_summary = pd.DataFrame(
        [
            {"metric": "final_species_filter", "value": receptor_result["final_species_filter"]},
            {"metric": "human_records_exist", "value": receptor_result["human_exists"]},
            {"metric": "distance_comparison_generated", "value": distance_result["generated"]},
        ]
    )

    backup_existing(workbook_path)
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        pairwise.to_excel(writer, sheet_name="receptor_pairwise_concordance", index=False)
        grouped.to_excel(writer, sheet_name="receptor_concordant_div_pairs", index=False)
        vectors.to_excel(writer, sheet_name="receptor_response_vectors", index=False)
        stats.to_excel(writer, sheet_name="receptor_consistency_statistics", index=False)
        distance_sheet.to_excel(writer, sheet_name="distance_comparison", index=False)
        audit_summary.to_excel(writer, sheet_name="audit_summary", index=False)
    return workbook_path


def main() -> int:
    receptor_result = prepare_m2or_receptor_annotations(PROJECT_ROOT)
    distance_result = compute_receptor_descriptor_distances(PROJECT_ROOT)
    _upsert_yaml_files(PROJECT_ROOT)
    _update_docs(PROJECT_ROOT)
    workbook_path = _build_fig4_workbook(PROJECT_ROOT, receptor_result, distance_result)

    print(f"receptor raw rows: {receptor_result['stats_map']['n_raw_rows']}")
    print(f"receptor pair count: {receptor_result['stats_map']['n_pairs_raw']}")
    print(f"pairs with jointly tested receptors: {receptor_result['stats_map']['n_pairs_with_jointly_tested_receptors']}")
    print(f"pairs with >=2 jointly tested receptors: {receptor_result['stats_map']['n_pairs_with_at_least_two_jointly_tested_receptors']}")
    print(f"concordant pairs: {receptor_result['stats_map']['n_concordant_pairs']}")
    print(f"divergent pairs: {receptor_result['stats_map']['n_divergent_pairs']}")
    print(f"indeterminate pairs: {receptor_result['stats_map']['n_indeterminate_pairs']}")
    print(f"excluded pairs: {receptor_result['stats_map']['n_excluded_pairs']}")
    print(f"descriptor-distance comparison generated: {distance_result['generated']}")
    print(f"fig4 source data: {workbook_path.relative_to(PROJECT_ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
