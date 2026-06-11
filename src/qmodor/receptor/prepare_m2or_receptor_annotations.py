"""Prepare M2OR receptor-response annotation tables from pairwise comparison data."""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from qmodor.descriptors.common import backup_existing, write_csv


ALLOWED_GROUPS = {
    "concordant",
    "divergent",
    "indeterminate",
    "insufficient_joint_receptor_data",
    "no_matched_or_data",
}


def snake_case(name: str) -> str:
    value = str(name).strip().lower()
    for char in [" ", "-", "/", "(", ")", ".", "%"]:
        value = value.replace(char, "_")
    while "__" in value:
        value = value.replace("__", "_")
    return value.strip("_")


def normalize_bool_response(value: Any) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in {"", "nan", "none", "n.d", "nd", "n/a", "not determined", "untested"}:
        return np.nan
    if text in {"1", "1.0", "true", "responsive", "yes", "active"}:
        return 1.0
    if text in {"0", "0.0", "false", "nonresponsive", "non-responsive", "no", "inactive"}:
        return 0.0
    try:
        number = float(text)
        if number == 1.0:
            return 1.0
        if number == 0.0:
            return 0.0
    except Exception:
        return np.nan
    return np.nan


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_species(value: Any) -> str:
    text = clean_text(value)
    if text.lower() in {"homo sapiens", "human"}:
        return "Homo sapiens"
    return text


def evidence_rank(value: Any) -> int:
    text = clean_text(value).lower()
    if text == "high":
        return 3
    if text == "medium":
        return 2
    if text == "low":
        return 1
    return 0


def _determine_pair_ids(values: pd.Series) -> dict[str, str]:
    unique_values = sorted({clean_text(v) for v in values if clean_text(v)}, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    return {value: f"PAIR_{idx:04d}" for idx, value in enumerate(unique_values, start=1)}


def _pair_group_label(n_joint: int, n_diff: int) -> tuple[str, bool]:
    if n_joint < 2:
        return "insufficient_joint_receptor_data", False
    if n_diff == 0:
        return "concordant", True
    if (n_diff / n_joint) >= 0.5:
        return "divergent", True
    return "indeterminate", False


def _upsert_markdown(path: Path, heading: str, body: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    section = f"\n## {heading}\n\n{body.strip()}\n"
    import re

    pattern = re.compile(rf"\n## {re.escape(heading)}\n.*?(?=\n## |\Z)", re.S)
    if existing:
        if pattern.search(existing):
            new_text = pattern.sub(section, existing)
        else:
            new_text = existing.rstrip() + section
    else:
        new_text = f"# {path.stem.replace('_', ' ').title()}\n" + section
    backup_existing(path)
    path.write_text(new_text.rstrip() + "\n", encoding="utf-8")


def prepare_m2or_receptor_annotations(project_root: Path) -> dict[str, Any]:
    """Prepare receptor annotation tables from the local M2OR pairwise comparison file."""

    workspace_root = project_root.parent
    input_path = workspace_root / "data_pairs_60_pair_OR_comparison.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Missing input file: {input_path}")

    receptor_dir = project_root / "data" / "receptor_annotations"
    results_dir = project_root / "results" / "receptor_consistency"
    reports_dir = project_root / "reports"
    source_data_dir = project_root / "data" / "source_data"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    source_data_dir.mkdir(parents=True, exist_ok=True)

    raw_input_df = pd.read_csv(input_path)
    raw_copy_output = receptor_dir / "m2or_pair_or_raw.csv"
    write_csv(raw_input_df, raw_copy_output)

    df = raw_input_df.copy()
    df.columns = [snake_case(col) for col in df.columns]
    if "pair_id_std" not in df.columns:
        raise ValueError("Input receptor comparison file is missing `pair_id_std`.")

    pair_mapping = _determine_pair_ids(df["pair_id_std"])
    df["pair_id_source"] = df["pair_id_std"].map(clean_text)
    df["pair_id"] = df["pair_id_source"].map(pair_mapping)
    df["cas_a"] = df["molecule_a_cas_norm"].map(clean_text)
    df.loc[df["cas_a"].eq(""), "cas_a"] = df.loc[df["cas_a"].eq(""), "molecule_a_cas_raw"].map(clean_text)
    df["cas_b"] = df["molecule_b_cas_norm"].map(clean_text)
    df.loc[df["cas_b"].eq(""), "cas_b"] = df.loc[df["cas_b"].eq(""), "molecule_b_cas_raw"].map(clean_text)
    df["pair_key"] = df.apply(lambda row: "__".join(sorted([row["cas_a"], row["cas_b"]])), axis=1)
    df["species"] = df["species"].map(normalize_species)
    df["or_gene"] = df["or_gene"].map(clean_text)
    df["a_responsive_norm"] = df["a_responsive"].map(normalize_bool_response)
    df["b_responsive_norm"] = df["b_responsive"].map(normalize_bool_response)
    df["a_tested"] = df["a_responsive_norm"].notna()
    df["b_tested"] = df["b_responsive_norm"].notna()
    df["jointly_tested"] = df["a_tested"] & df["b_tested"]
    df["same_response"] = np.where(
        df["jointly_tested"],
        df["a_responsive_norm"] == df["b_responsive_norm"],
        np.nan,
    )
    df["differential_response"] = np.where(
        df["jointly_tested"],
        df["a_responsive_norm"] != df["b_responsive_norm"],
        np.nan,
    )
    df["base_exclusion_reason"] = ""
    df.loc[df["or_gene"].eq(""), "base_exclusion_reason"] = "missing_or_gene"
    df.loc[df["or_gene"].eq("") & df["comparison_type"].fillna("").astype(str).str.contains("No matched OR data", case=False), "base_exclusion_reason"] = "no_matched_or_data"
    df.loc[df["base_exclusion_reason"].eq("") & ~df["a_tested"] & ~df["b_tested"], "base_exclusion_reason"] = "both_untested"
    df.loc[df["base_exclusion_reason"].eq("") & ~df["a_tested"] & df["b_tested"], "base_exclusion_reason"] = "a_untested"
    df.loc[df["base_exclusion_reason"].eq("") & df["a_tested"] & ~df["b_tested"], "base_exclusion_reason"] = "b_untested"
    df["evidence_rank"] = df["evidence_priority"].map(evidence_rank)
    df["doi_presence_score"] = (
        df["a_doi"].fillna("").astype(str).str.strip().ne("").astype(int)
        + df["b_doi"].fillna("").astype(str).str.strip().ne("").astype(int)
    )

    pair_meta = (
        df[
            [
                "pair_id",
                "pair_id_source",
                "pair_key",
                "cas_a",
                "cas_b",
                "molecule_a_name",
                "molecule_b_name",
                "molecule_a_smiles",
                "molecule_b_smiles",
                "is_must_keep",
            ]
        ]
        .drop_duplicates(subset=["pair_id"])
        .sort_values("pair_id")
        .reset_index(drop=True)
    )

    harmonized_rows = []
    exclusion_rows = []
    duplicate_conflicts = 0
    duplicate_groups = df.groupby(["pair_id", "species", "or_gene"], dropna=False)
    for (_, species, or_gene), group in duplicate_groups:
        work = group.copy()
        if clean_text(or_gene) == "":
            for _, row in work.iterrows():
                exclusion_rows.append(
                    {
                        "pair_id": row["pair_id"],
                        "pair_key": row["pair_key"],
                        "species": row["species"],
                        "or_gene": row["or_gene"],
                        "row_status": "excluded",
                        "exclusion_reason": row["base_exclusion_reason"] or "missing_or_gene",
                        "pair_id_source": row["pair_id_source"],
                        "notes": row.get("notes", ""),
                    }
                )
            continue

        work = work.sort_values(
            by=["jointly_tested", "evidence_rank", "doi_presence_score"],
            ascending=[False, False, False],
            kind="mergesort",
        ).reset_index(drop=True)
        top = work.iloc[0]
        top_mask = (
            (work["jointly_tested"] == top["jointly_tested"])
            & (work["evidence_rank"] == top["evidence_rank"])
            & (work["doi_presence_score"] == top["doi_presence_score"])
        )
        top_group = work.loc[top_mask].copy()
        signatures = {
            (
                "nan" if pd.isna(row["a_responsive_norm"]) else int(row["a_responsive_norm"]),
                "nan" if pd.isna(row["b_responsive_norm"]) else int(row["b_responsive_norm"]),
            )
            for _, row in top_group.iterrows()
        }
        if len(signatures) > 1:
            duplicate_conflicts += 1
            for _, row in work.iterrows():
                exclusion_rows.append(
                    {
                        "pair_id": row["pair_id"],
                        "pair_key": row["pair_key"],
                        "species": row["species"],
                        "or_gene": row["or_gene"],
                        "row_status": "excluded",
                        "exclusion_reason": "conflicting_duplicate_records",
                        "pair_id_source": row["pair_id_source"],
                        "notes": row.get("notes", ""),
                    }
                )
            continue

        selected = top_group.iloc[0].copy()
        harmonized_rows.append(selected)
        for _, row in work.iloc[1:].iterrows():
            exclusion_rows.append(
                {
                    "pair_id": row["pair_id"],
                    "pair_key": row["pair_key"],
                    "species": row["species"],
                    "or_gene": row["or_gene"],
                    "row_status": "excluded",
                    "exclusion_reason": "lower_priority_duplicate_record",
                    "pair_id_source": row["pair_id_source"],
                    "notes": row.get("notes", ""),
                }
            )

    harmonized_df = pd.DataFrame(harmonized_rows)
    if harmonized_df.empty:
        harmonized_df = pd.DataFrame(columns=df.columns.tolist())

    human_exists = bool((harmonized_df["species"] == "Homo sapiens").any()) if not harmonized_df.empty else False
    final_species_filter = "human_only" if human_exists else "all_species"
    final_df = harmonized_df.loc[harmonized_df["species"] == "Homo sapiens"].copy() if human_exists else harmonized_df.copy()

    if human_exists:
        nonhuman = harmonized_df.loc[harmonized_df["species"] != "Homo sapiens"].copy()
        for _, row in nonhuman.iterrows():
            exclusion_rows.append(
                {
                    "pair_id": row["pair_id"],
                    "pair_key": row["pair_key"],
                    "species": row["species"],
                    "or_gene": row["or_gene"],
                    "row_status": "excluded",
                    "exclusion_reason": "non_human_species_if_human_filter_applied",
                    "pair_id_source": row["pair_id_source"],
                    "notes": row.get("notes", ""),
                }
            )

    for frame in [harmonized_df, final_df]:
        if frame.empty:
            continue
        frame["same_response"] = frame["same_response"].astype("float")
        frame["differential_response"] = frame["differential_response"].astype("float")

    def build_long_table(frame: pd.DataFrame, species_filter_label: str) -> pd.DataFrame:
        records = []
        for _, row in frame.iterrows():
            exclusion_reason = row["base_exclusion_reason"]
            records.append(
                {
                    "pair_id": row["pair_id"],
                    "pair_key": row["pair_key"],
                    "cas_a": row["cas_a"],
                    "cas_b": row["cas_b"],
                    "smiles_a": row["molecule_a_smiles"],
                    "smiles_b": row["molecule_b_smiles"],
                    "molecule_a_name": row["molecule_a_name"],
                    "molecule_b_name": row["molecule_b_name"],
                    "species": row["species"],
                    "or_gene": row["or_gene"],
                    "main_receptors_id": row.get("main_receptors_id", np.nan),
                    "a_responsive": row["a_responsive_norm"],
                    "b_responsive": row["b_responsive_norm"],
                    "a_tested": bool(row["a_tested"]),
                    "b_tested": bool(row["b_tested"]),
                    "jointly_tested": bool(row["jointly_tested"]),
                    "same_response": row["same_response"] if bool(row["jointly_tested"]) else np.nan,
                    "differential_response": row["differential_response"] if bool(row["jointly_tested"]) else np.nan,
                    "comparison_type": row.get("comparison_type", np.nan),
                    "evidence_priority": row.get("evidence_priority", np.nan),
                    "a_doi": row.get("a_doi", np.nan),
                    "b_doi": row.get("b_doi", np.nan),
                    "notes": row.get("notes", np.nan),
                    "exclusion_reason": exclusion_reason,
                    "species_filter": species_filter_label,
                    "pair_id_source": row["pair_id_source"],
                }
            )
        return pd.DataFrame(records)

    long_all_species = build_long_table(harmonized_df, "all_species")
    long_final = build_long_table(final_df, final_species_filter)

    def build_pairwise_concordance(frame: pd.DataFrame, species_filter_label: str) -> pd.DataFrame:
        rows = []
        frame_by_pair = frame.groupby("pair_id", dropna=False) if not frame.empty else {}
        for _, meta_row in pair_meta.iterrows():
            pair_id = meta_row["pair_id"]
            if not frame.empty and pair_id in frame_by_pair.groups:
                pair_frame = frame_by_pair.get_group(pair_id).copy()
                n_or_rows_raw = int(len(pair_frame))
                jointly = pair_frame.loc[pair_frame["jointly_tested"]]
                n_joint = int(len(jointly))
                n_same = int((jointly["same_response"] == 1).sum())
                n_diff = int((jointly["differential_response"] == 1).sum())
                same_prop = (n_same / n_joint) if n_joint else np.nan
                diff_prop = (n_diff / n_joint) if n_joint else np.nan
                group_label, include = _pair_group_label(n_joint, n_diff)
                exclusion_reason = ""
                if n_joint < 2:
                    exclusion_reason = "insufficient_jointly_tested_receptors"
                rows.append(
                    {
                        "pair_id": pair_id,
                        "pair_key": meta_row["pair_key"],
                        "cas_a": meta_row["cas_a"],
                        "cas_b": meta_row["cas_b"],
                        "n_or_rows_raw": n_or_rows_raw,
                        "n_jointly_tested_receptors": n_joint,
                        "n_same_response": n_same,
                        "n_differential_response": n_diff,
                        "same_response_proportion": same_prop,
                        "differential_response_proportion": diff_prop,
                        "receptor_response_group": group_label,
                        "included_in_group_comparison": bool(include),
                        "exclusion_reason": exclusion_reason,
                        "jointly_tested_or_genes": ";".join(sorted(jointly["or_gene"].dropna().astype(str).unique().tolist())),
                        "same_response_or_genes": ";".join(sorted(jointly.loc[jointly["same_response"] == 1, "or_gene"].dropna().astype(str).unique().tolist())),
                        "differential_response_or_genes": ";".join(sorted(jointly.loc[jointly["differential_response"] == 1, "or_gene"].dropna().astype(str).unique().tolist())),
                        "evidence_summary": ";".join(sorted(pair_frame["evidence_priority"].fillna("").astype(str).replace("", pd.NA).dropna().unique().tolist())),
                        "species_filter": species_filter_label,
                    }
                )
            else:
                exclusion_reason = "no_matched_or_data"
                if human_exists and not harmonized_df.empty and pair_id in set(harmonized_df["pair_id"]):
                    exclusion_reason = "non_human_species_if_human_filter_applied"
                rows.append(
                    {
                        "pair_id": pair_id,
                        "pair_key": meta_row["pair_key"],
                        "cas_a": meta_row["cas_a"],
                        "cas_b": meta_row["cas_b"],
                        "n_or_rows_raw": 0,
                        "n_jointly_tested_receptors": 0,
                        "n_same_response": 0,
                        "n_differential_response": 0,
                        "same_response_proportion": np.nan,
                        "differential_response_proportion": np.nan,
                        "receptor_response_group": "no_matched_or_data",
                        "included_in_group_comparison": False,
                        "exclusion_reason": exclusion_reason,
                        "jointly_tested_or_genes": "",
                        "same_response_or_genes": "",
                        "differential_response_or_genes": "",
                        "evidence_summary": "",
                        "species_filter": species_filter_label,
                    }
                )
        return pd.DataFrame(rows)

    pairwise_all_species = build_pairwise_concordance(harmonized_df, "all_species")
    pairwise_final = build_pairwise_concordance(final_df, final_species_filter)

    included_final = pairwise_final.loc[pairwise_final["included_in_group_comparison"]].copy()
    included_final = included_final.loc[included_final["receptor_response_group"].isin(["concordant", "divergent"])].reset_index(drop=True)

    def build_matched_molecules(frame: pd.DataFrame, species_filter_label: str) -> pd.DataFrame:
        molecule_rows = []
        for _, row in frame.iterrows():
            molecule_rows.append(
                {
                    "cas": row["cas_a"],
                    "molecule_role_source": "A",
                    "molecule_name": row["molecule_a_name"],
                    "smiles": row["molecule_a_smiles"],
                    "m2or_main_compounds_id": row.get("molecule_a_main_compounds_id", np.nan),
                    "pair_id": row["pair_id"],
                    "or_gene": row["or_gene"],
                    "responsive": row["a_responsive_norm"],
                    "jointly_tested": bool(row["jointly_tested"]),
                    "species_filter": species_filter_label,
                }
            )
            molecule_rows.append(
                {
                    "cas": row["cas_b"],
                    "molecule_role_source": "B",
                    "molecule_name": row["molecule_b_name"],
                    "smiles": row["molecule_b_smiles"],
                    "m2or_main_compounds_id": row.get("molecule_b_main_compounds_id", np.nan),
                    "pair_id": row["pair_id"],
                    "or_gene": row["or_gene"],
                    "responsive": row["b_responsive_norm"],
                    "jointly_tested": bool(row["jointly_tested"]),
                    "species_filter": species_filter_label,
                }
            )
        molecules_df = pd.DataFrame(molecule_rows)
        if molecules_df.empty:
            return pd.DataFrame(
                columns=[
                    "cas",
                    "molecule_role_source",
                    "molecule_name",
                    "smiles",
                    "m2or_main_compounds_id",
                    "n_pairs",
                    "n_receptor_rows",
                    "n_jointly_tested_rows",
                    "n_responsive_records",
                    "n_nonresponsive_records",
                ]
            )
        summary = (
            molecules_df.groupby(["cas", "molecule_role_source", "molecule_name", "smiles", "m2or_main_compounds_id"], dropna=False)
            .agg(
                n_pairs=("pair_id", "nunique"),
                n_receptor_rows=("or_gene", "count"),
                n_jointly_tested_rows=("jointly_tested", "sum"),
                n_responsive_records=("responsive", lambda s: int((s == 1).sum())),
                n_nonresponsive_records=("responsive", lambda s: int((s == 0).sum())),
            )
            .reset_index()
        )
        return summary

    matched_all_species = build_matched_molecules(harmonized_df, "all_species")
    matched_final = build_matched_molecules(final_df, final_species_filter)

    def build_response_vectors(frame: pd.DataFrame, species_filter_label: str) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
        records = []
        vector_exclusions = []
        for _, row in frame.iterrows():
            for role, cas_col, name_col, smiles_col, response_col in [
                ("A", "cas_a", "molecule_a_name", "molecule_a_smiles", "a_responsive_norm"),
                ("B", "cas_b", "molecule_b_name", "molecule_b_smiles", "b_responsive_norm"),
            ]:
                records.append(
                    {
                        "cas": row[cas_col],
                        "molecule_name": row[name_col],
                        "smiles": row[smiles_col],
                        "species_filter": species_filter_label,
                        "or_gene": row["or_gene"],
                        "responsive": row[response_col],
                        "pair_id": row["pair_id"],
                    }
                )
        molecule_or_df = pd.DataFrame(records)
        if molecule_or_df.empty:
            return pd.DataFrame(columns=["cas", "molecule_name", "smiles", "species_filter"]), vector_exclusions

        agg_rows = []
        for (cas, molecule_name, smiles), group in molecule_or_df.groupby(["cas", "molecule_name", "smiles"], dropna=False):
            row = {"cas": cas, "molecule_name": molecule_name, "smiles": smiles, "species_filter": species_filter_label}
            for or_gene, receptor_group in group.groupby("or_gene", dropna=False):
                values = receptor_group["responsive"].dropna().unique().tolist()
                if len(values) == 1:
                    row[or_gene] = int(values[0])
                elif len(values) == 0:
                    row[or_gene] = np.nan
                else:
                    row[or_gene] = np.nan
                    vector_exclusions.append(
                        {
                            "pair_id": ";".join(sorted(receptor_group["pair_id"].astype(str).unique().tolist())),
                            "pair_key": "",
                            "species": species_filter_label,
                            "or_gene": or_gene,
                            "row_status": "excluded",
                            "exclusion_reason": "conflicting_duplicate_records",
                            "pair_id_source": "",
                            "notes": f"Conflicting molecule-level response vector values for CAS {cas}.",
                        }
                    )
            agg_rows.append(row)
        return pd.DataFrame(agg_rows), vector_exclusions

    vectors_all_species, vector_exclusions_all = build_response_vectors(harmonized_df, "all_species")
    vectors_final, vector_exclusions_final = build_response_vectors(final_df, final_species_filter)
    exclusion_rows.extend(vector_exclusions_all)
    exclusion_rows.extend(vector_exclusions_final)

    exclusion_df = pd.DataFrame(exclusion_rows).drop_duplicates().reset_index(drop=True)

    stats_rows = []
    species_distribution = df["species"].replace("", pd.NA).fillna("missing").value_counts(dropna=False).to_dict()
    evidence_distribution = df["evidence_priority"].replace("", pd.NA).fillna("missing").value_counts(dropna=False).to_dict()
    jointly_distribution = pairwise_final["n_jointly_tested_receptors"].value_counts(dropna=False).sort_index().to_dict()
    stats_map = {
        "n_raw_rows": int(len(df)),
        "n_pairs_raw": int(pair_meta["pair_id"].nunique()),
        "n_pairs_with_any_or_record": int((pairwise_final["n_or_rows_raw"] > 0).sum()),
        "n_pairs_no_matched_or_data": int((pairwise_final["receptor_response_group"] == "no_matched_or_data").sum()),
        "n_pairs_with_jointly_tested_receptors": int((pairwise_final["n_jointly_tested_receptors"] > 0).sum()),
        "n_pairs_with_at_least_two_jointly_tested_receptors": int((pairwise_final["n_jointly_tested_receptors"] >= 2).sum()),
        "n_concordant_pairs": int((pairwise_final["receptor_response_group"] == "concordant").sum()),
        "n_divergent_pairs": int((pairwise_final["receptor_response_group"] == "divergent").sum()),
        "n_indeterminate_pairs": int((pairwise_final["receptor_response_group"] == "indeterminate").sum()),
        "n_excluded_pairs": int((pairwise_final["included_in_group_comparison"] == False).sum()),
        "n_unique_or_genes": int(df["or_gene"].replace("", pd.NA).dropna().nunique()),
        "n_unique_human_or_genes": int(harmonized_df.loc[harmonized_df["species"] == "Homo sapiens", "or_gene"].replace("", pd.NA).dropna().nunique()),
        "n_conflicting_duplicate_records": int(duplicate_conflicts),
    }
    for metric, value in stats_map.items():
        stats_rows.append({"metric": metric, "value": value})
    stats_df = pd.DataFrame(stats_rows)

    outputs = {
        "m2or_pair_or_raw": raw_copy_output,
        "m2or_matched_molecules": receptor_dir / "m2or_matched_molecules.csv",
        "m2or_matched_molecules_all_species": receptor_dir / "m2or_matched_molecules_all_species.csv",
        "receptor_response_long": receptor_dir / "receptor_response_long.csv",
        "receptor_response_long_all_species": receptor_dir / "receptor_response_long_all_species.csv",
        "receptor_response_vectors": receptor_dir / "receptor_response_vectors.csv",
        "receptor_response_vectors_all_species": receptor_dir / "receptor_response_vectors_all_species.csv",
        "receptor_pairwise_concordance": receptor_dir / "receptor_pairwise_concordance.csv",
        "receptor_pairwise_concordance_all_species": receptor_dir / "receptor_pairwise_concordance_all_species.csv",
        "receptor_concordant_divergent_pairs": receptor_dir / "receptor_concordant_divergent_pairs.csv",
        "receptor_concordant_divergent_pairs_all_species": receptor_dir / "receptor_concordant_divergent_pairs_all_species.csv",
        "receptor_matching_exclusion_log": receptor_dir / "receptor_matching_exclusion_log.csv",
        "receptor_consistency_statistics": results_dir / "receptor_consistency_statistics.csv",
    }

    write_csv(matched_final, outputs["m2or_matched_molecules"])
    write_csv(matched_all_species, outputs["m2or_matched_molecules_all_species"])
    write_csv(long_final, outputs["receptor_response_long"])
    write_csv(long_all_species, outputs["receptor_response_long_all_species"])
    write_csv(vectors_final, outputs["receptor_response_vectors"])
    write_csv(vectors_all_species, outputs["receptor_response_vectors_all_species"])
    write_csv(pairwise_final, outputs["receptor_pairwise_concordance"])
    write_csv(pairwise_all_species, outputs["receptor_pairwise_concordance_all_species"])
    write_csv(included_final, outputs["receptor_concordant_divergent_pairs"])
    write_csv(
        pairwise_all_species.loc[
            pairwise_all_species["included_in_group_comparison"]
            & pairwise_all_species["receptor_response_group"].isin(["concordant", "divergent"])
        ].reset_index(drop=True),
        outputs["receptor_concordant_divergent_pairs_all_species"],
    )
    write_csv(exclusion_df, outputs["receptor_matching_exclusion_log"])
    write_csv(stats_df, outputs["receptor_consistency_statistics"])

    report_lines = [
        "# Receptor Annotation Preparation Report",
        "",
        "## 1. Input file inspected",
        f"- `{input_path.name}`",
        "",
        "## 2. Raw row count",
        f"- `{len(df)}`",
        "",
        "## 3. Pair count",
        f"- `{pair_meta['pair_id'].nunique()}`",
        "",
        "## 4. Unique CAS count",
        f"- `{len(set(df['cas_a']) | set(df['cas_b']))}`",
        "",
        "## 5. Unique OR_gene count",
        f"- `{stats_map['n_unique_or_genes']}`",
        "",
        "## 6. Species distribution",
    ]
    report_lines.extend([f"- `{species}`: `{count}`" for species, count in species_distribution.items()])
    report_lines.extend(["", "## 7. Evidence priority distribution"])
    report_lines.extend([f"- `{level}`: `{count}`" for level, count in evidence_distribution.items()])
    report_lines.extend(["", "## 8. Jointly tested receptor count distribution by pair"])
    report_lines.extend([f"- `{key}` jointly tested receptors: `{value}` pairs" for key, value in jointly_distribution.items()])
    report_lines.extend(
        [
            "",
            "## 9. Number of concordant/divergent/indeterminate/excluded pairs",
            f"- concordant: `{stats_map['n_concordant_pairs']}`",
            f"- divergent: `{stats_map['n_divergent_pairs']}`",
            f"- indeterminate: `{stats_map['n_indeterminate_pairs']}`",
            f"- excluded: `{stats_map['n_excluded_pairs']}`",
            "",
            "## 10. Duplicate/conflict handling",
            f"- conflicting duplicate record groups: `{duplicate_conflicts}`",
            "- Evidence priority ranking used: High > Medium > Low > missing.",
            "- Jointly tested records were preferred over untested records.",
            "- DOI-bearing records were preferred over records without DOI when other ranks tied.",
            "",
            "## 11. Files generated",
        ]
    )
    report_lines.extend([f"- `{path.relative_to(project_root).as_posix()}`" for path in outputs.values()])
    report_lines.extend(["", "## 12. Files not generated and why", "- Descriptor-distance comparison files are handled separately and may be blocked if descriptor matrices are unavailable.", "", "## 13. Interpretation boundary", "These public receptor-response annotations are used as an orthogonal consistency analysis and do not establish receptor-resolved binding mechanisms."])
    report_path = reports_dir / "receptor_annotation_preparation_report.md"
    from qmodor.descriptors.common import write_text

    write_text("\n".join(report_lines) + "\n", report_path)

    receptor_rules_path = project_root / "docs" / "receptor_matching_rules.md"
    receptor_rules_body = "\n".join(
        [
            f"- Current local input file: `{input_path.name}`.",
            "- Column names are normalized to snake_case in processed outputs.",
            "- Pair IDs are standardized to `PAIR_0001` style while preserving the original `pair_id_std` in `pair_id_source`.",
            "- CAS matching defaults to normalized CAS columns and falls back to raw CAS when needed.",
            "- Tested status is defined by non-missing responsive annotation after cautious normalization of 1/0/True/False/responsive/nonresponsive strings.",
            "- Duplicate pair-by-species-by-OR records are ranked by jointly tested status, evidence priority, and DOI presence; unresolved conflicts are logged and excluded.",
            "- Concordant pairs have zero differential responses across jointly tested receptors; divergent pairs have differential-response proportion greater than or equal to 0.5; indeterminate pairs are excluded from grouped comparison; fewer than two jointly tested receptors are marked insufficient.",
            "- Human-only outputs are used as the default manuscript-consistent view when Homo sapiens records exist; all-species companion outputs are also retained.",
            "- These public receptor-response annotations are used as an orthogonal consistency analysis and do not establish receptor-resolved binding mechanisms.",
        ]
    )
    _upsert_markdown(receptor_rules_path, "Current Local Processing Rules", receptor_rules_body)

    return {
        "raw_df": df,
        "pair_meta": pair_meta,
        "long_final": long_final,
        "pairwise_final": pairwise_final,
        "included_final": included_final,
        "stats_df": stats_df,
        "stats_map": stats_map,
        "species_distribution": species_distribution,
        "evidence_distribution": evidence_distribution,
        "human_exists": human_exists,
        "outputs": outputs,
        "report_path": report_path,
        "final_species_filter": final_species_filter,
    }
