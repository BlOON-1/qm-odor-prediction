"""Compute raw RDKit 2D descriptors from the benchmark molecule table."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import (
    InputTableInfo,
    clean_numeric_frame,
    load_benchmark_table,
    require_rdkit,
    resolve_smiles_column,
    write_csv,
)


LOGGER = logging.getLogger(__name__)


def _canonicalize_or_none(chem: Any, smiles: str) -> tuple[Any, str]:
    """Create an RDKit molecule and canonical SMILES when possible."""

    if not smiles:
        return None, ""
    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        return None, ""
    return mol, chem.MolToSmiles(mol, canonical=True)


def compute_rdkit_descriptors(
    project_root: Path,
    input_path: str | None = None,
) -> dict[str, Any]:
    """Compute RDKit descriptors and audit outputs for valid benchmark molecules."""

    chem, descriptor_module, rdkit_pkg = require_rdkit()
    input_info: InputTableInfo = load_benchmark_table(project_root, input_path)
    benchmark_df = resolve_smiles_column(input_info.df)

    reports_dir = project_root / "reports"
    processed_dir = project_root / "data" / "processed"
    reports_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    descriptor_names = [name for name, _ in descriptor_module._descList]
    descriptor_fns = dict(descriptor_module._descList)

    invalid_rows: list[dict[str, Any]] = []
    valid_rows: list[dict[str, Any]] = []
    descriptor_failure_counts = {name: 0 for name in descriptor_names}

    for _, row in benchmark_df.iterrows():
        source_smiles = str(row["smiles_input"]).strip()
        mol, canonical_smiles = _canonicalize_or_none(chem, source_smiles)
        if mol is None:
            invalid_rows.append(
                {
                    "molecule_id": row["molecule_id"],
                    "cas": row["cas"],
                    "chemical_name": row["chemical_name"],
                    "smiles_input": source_smiles,
                    "failure_reason": "invalid_or_missing_smiles",
                }
            )
            continue

        record = {
            "molecule_id": row["molecule_id"],
            "cas": row["cas"],
            "chemical_name": row["chemical_name"],
            "canonical_smiles": canonical_smiles,
        }
        for descriptor_name in descriptor_names:
            descriptor_fn = descriptor_fns[descriptor_name]
            try:
                record[descriptor_name] = descriptor_fn(mol)
            except Exception:
                record[descriptor_name] = np.nan
                descriptor_failure_counts[descriptor_name] += 1
        valid_rows.append(record)

    raw_df = pd.DataFrame(valid_rows)
    invalid_df = pd.DataFrame(
        invalid_rows,
        columns=["molecule_id", "cas", "chemical_name", "smiles_input", "failure_reason"],
    )
    if raw_df.empty:
        raw_df = pd.DataFrame(columns=["molecule_id", "cas", "chemical_name", "canonical_smiles", *descriptor_names])

    metadata_cols = ["molecule_id", "cas", "chemical_name", "canonical_smiles"]
    numeric_cols = [col for col in raw_df.columns if col not in metadata_cols]
    cleaned_numeric, inf_counts = clean_numeric_frame(raw_df[numeric_cols]) if numeric_cols else (pd.DataFrame(), pd.Series(dtype="int64"))
    if numeric_cols:
        raw_df.loc[:, numeric_cols] = cleaned_numeric

    missingness_rows = []
    all_missing_cols = []
    constant_rows = []
    for descriptor_name in numeric_cols:
        column = raw_df[descriptor_name]
        missing_count = int(column.isna().sum())
        missing_fraction = float(missing_count / len(raw_df)) if len(raw_df) else np.nan
        n_unique_non_na = int(column.dropna().nunique())
        if missing_count == len(raw_df):
            all_missing_cols.append(descriptor_name)
        if n_unique_non_na <= 1 and missing_count < len(raw_df):
            constant_rows.append(
                {
                    "descriptor": descriptor_name,
                    "constant_value": column.dropna().iloc[0] if n_unique_non_na == 1 else np.nan,
                    "n_unique_non_na": n_unique_non_na,
                }
            )
        missingness_rows.append(
            {
                "descriptor": descriptor_name,
                "missing_count": missing_count,
                "missing_fraction": missing_fraction,
                "inf_replaced_count": int(inf_counts.get(descriptor_name, 0)),
                "calculation_failure_count": int(descriptor_failure_counts.get(descriptor_name, 0)),
            }
        )

    filtered_raw_df = raw_df.drop(columns=all_missing_cols, errors="ignore")
    missingness_df = pd.DataFrame(missingness_rows).sort_values("descriptor").reset_index(drop=True)
    constant_df = pd.DataFrame(
        constant_rows,
        columns=["descriptor", "constant_value", "n_unique_non_na"],
    ).sort_values("descriptor").reset_index(drop=True)
    audit_df = pd.DataFrame(
        [
            {
                "input_file": str(input_info.path.relative_to(project_root)),
                "rdkit_version": rdkit_pkg.__version__,
                "n_molecules_total": int(len(benchmark_df)),
                "n_molecules_valid": int(len(raw_df)),
                "n_molecules_invalid": int(len(invalid_df)),
                "n_descriptor_columns_before_all_missing_removal": int(len(numeric_cols)),
                "n_descriptor_columns_after_all_missing_removal": int(len(filtered_raw_df.columns) - len(metadata_cols)),
                "n_all_missing_descriptor_columns_removed": int(len(all_missing_cols)),
                "total_inf_values_replaced_with_nan": int(inf_counts.sum()) if len(inf_counts) else 0,
            }
        ]
    )

    raw_output = processed_dir / "rdkit_descriptors_raw.csv"
    invalid_output = reports_dir / "rdkit_invalid_smiles.csv"
    audit_output = reports_dir / "rdkit_descriptor_audit.csv"
    missingness_output = reports_dir / "rdkit_descriptor_missingness.csv"
    constant_output = reports_dir / "rdkit_descriptor_constant_columns.csv"

    write_csv(filtered_raw_df, raw_output)
    write_csv(invalid_df, invalid_output)
    write_csv(audit_df, audit_output)
    write_csv(missingness_df, missingness_output)
    write_csv(constant_df, constant_output)

    LOGGER.info(
        "Computed RDKit descriptors for %s valid molecules; %s invalid SMILES entries were logged.",
        len(raw_df),
        len(invalid_df),
    )
    return {
        "input_info": input_info,
        "raw_df": filtered_raw_df,
        "invalid_df": invalid_df,
        "audit_df": audit_df,
        "missingness_df": missingness_df,
        "constant_df": constant_df,
        "descriptor_names_before_filtering": descriptor_names,
        "all_missing_removed": all_missing_cols,
        "output_files": [
            raw_output,
            invalid_output,
            audit_output,
            missingness_output,
            constant_output,
        ],
    }
