"""Train-split-aware standardization for RDKit continuous descriptors."""

from __future__ import annotations

import itertools
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import require_sklearn, write_csv, write_text


LOGGER = logging.getLogger(__name__)

METADATA_COLUMNS = ["molecule_id", "cas", "chemical_name", "canonical_smiles"]


def _discover_split_pairs(project_root: Path) -> list[tuple[int, Path, Path]]:
    """Locate seed-specific train/test split files."""

    split_dir = project_root / "data" / "splits" / "iterative_stratified_splits"
    pairs: list[tuple[int, Path, Path]] = []
    for seed in range(42, 47):
        train_path = split_dir / f"seed_{seed}_train.csv"
        test_path = split_dir / f"seed_{seed}_test.csv"
        if train_path.exists() and test_path.exists():
            pairs.append((seed, train_path, test_path))
    return pairs


def _load_split_ids(path: Path) -> pd.Series:
    """Read molecule IDs from a split file."""

    df = pd.read_csv(path)
    if "molecule_id" not in df.columns:
        raise ValueError(f"Split file {path} does not contain `molecule_id`.")
    return df["molecule_id"].astype(str)


def _rank_features_for_correlation(train_df: pd.DataFrame, missingness: pd.Series) -> list[str]:
    """Rank features so stronger candidates are kept during correlation pruning."""

    variance = train_df.var(axis=0, ddof=0)
    ranked = sorted(
        train_df.columns.tolist(),
        key=lambda col: (
            float(missingness.get(col, 0.0)),
            -float(variance.get(col, 0.0)),
            col,
        ),
    )
    return ranked


def _apply_training_filters(
    train_features: pd.DataFrame,
    test_features: pd.DataFrame,
    threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """Apply missingness, constant-column, and correlation filtering using training data only."""

    filtering_log_rows = []
    missingness = train_features.isna().mean(axis=0)
    missing_columns = sorted(missingness[missingness > 0].index.tolist())
    for column in missing_columns:
        filtering_log_rows.append(
            {
                "feature": column,
                "action": "removed",
                "stage": "missing_value_filter",
                "reason": "contains_missing_values_in_training_set",
                "metric_value": float(missingness[column]),
            }
        )

    train_step = train_features.drop(columns=missing_columns, errors="ignore")
    test_step = test_features.drop(columns=missing_columns, errors="ignore")

    variance = train_step.var(axis=0, ddof=0)
    constant_columns = sorted(variance[variance == 0].index.tolist())
    for column in constant_columns:
        filtering_log_rows.append(
            {
                "feature": column,
                "action": "removed",
                "stage": "constant_filter",
                "reason": "zero_variance_in_training_set",
                "metric_value": 0.0,
            }
        )
    train_step = train_step.drop(columns=constant_columns, errors="ignore")
    test_step = test_step.drop(columns=constant_columns, errors="ignore")

    retained = _rank_features_for_correlation(train_step, missingness)
    corr = train_step[retained].corr(method="pearson").abs() if retained else pd.DataFrame()
    kept: list[str] = []
    removed_due_to_corr: list[str] = []
    for feature in retained:
        correlated_with_kept = [
            keep_feature
            for keep_feature in kept
            if float(corr.loc[feature, keep_feature]) > threshold
        ]
        if correlated_with_kept:
            removed_due_to_corr.append(feature)
            max_corr = max(float(corr.loc[feature, keep_feature]) for keep_feature in correlated_with_kept)
            filtering_log_rows.append(
                {
                    "feature": feature,
                    "action": "removed",
                    "stage": "correlation_filter",
                    "reason": "high_absolute_pearson_correlation",
                    "metric_value": max_corr,
                    "correlated_with_retained_feature": correlated_with_kept[0],
                }
            )
        else:
            kept.append(feature)
            filtering_log_rows.append(
                {
                    "feature": feature,
                    "action": "retained",
                    "stage": "correlation_filter",
                    "reason": "passed_correlation_filter",
                    "metric_value": np.nan,
                }
            )

    train_final = train_step[kept].copy()
    test_final = test_step[kept].copy()
    filtering_log = pd.DataFrame(filtering_log_rows)
    return train_final, test_final, filtering_log, kept


def standardize_rdkit_descriptors(
    project_root: Path,
    correlation_threshold: float = 0.9,
    allow_global_qc_standardization: bool = False,
) -> dict[str, Any]:
    """Standardize RDKit descriptors using train-set parameters only when splits exist."""

    StandardScaler = require_sklearn()
    raw_path = project_root / "data" / "processed" / "rdkit_descriptors_raw.csv"
    if not raw_path.exists():
        raise FileNotFoundError(f"Expected raw RDKit descriptor file was not found: {raw_path}")

    raw_df = pd.read_csv(raw_path)
    metadata_cols = [col for col in METADATA_COLUMNS if col in raw_df.columns]
    feature_cols = [col for col in raw_df.columns if col not in metadata_cols]
    features = raw_df[feature_cols].copy()

    split_pairs = _discover_split_pairs(project_root)
    processed_dir = project_root / "data" / "processed"
    standardized_dir = processed_dir / "standardized"
    reports_dir = project_root / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    missing_columns_global = sorted(features.columns[features.isna().any(axis=0)].tolist())
    qc_df = raw_df.drop(columns=missing_columns_global, errors="ignore")
    constant_cols_global = []
    global_feature_cols = [col for col in qc_df.columns if col not in metadata_cols]
    if global_feature_cols:
        variance_global = qc_df[global_feature_cols].var(axis=0, ddof=0)
        constant_cols_global = sorted(variance_global[variance_global == 0].index.tolist())
        qc_df = qc_df.drop(columns=constant_cols_global, errors="ignore")

    qc_output = processed_dir / "rdkit_descriptors_unstandardized_qc.csv"
    write_csv(qc_df, qc_output)

    output_files = [qc_output]
    retained_feature_counts: dict[int, int] = {}
    generated_standardized = False

    if not split_pairs:
        blocked_report = (
            "# RDKit Standardization Blocked\n\n"
            "Final standardized RDKit descriptors were not generated because train/test split files were not available. "
            "Fitting a scaler on all molecules would introduce data leakage and is therefore not performed.\n"
        )
        blocked_path = reports_dir / "rdkit_standardization_blocked_no_splits.md"
        write_text(blocked_report, blocked_path)
        output_files.append(blocked_path)

        if allow_global_qc_standardization:
            scaler = StandardScaler()
            global_standardized = raw_df[metadata_cols].copy()
            global_standardized_features = scaler.fit_transform(qc_df[[col for col in qc_df.columns if col not in metadata_cols]])
            global_standardized = pd.concat(
                [
                    global_standardized,
                    pd.DataFrame(
                        global_standardized_features,
                        columns=[col for col in qc_df.columns if col not in metadata_cols],
                    ),
                ],
                axis=1,
            )
            global_qc_output = processed_dir / "rdkit_descriptors_global_standardized_FOR_QC_ONLY_NOT_FOR_MODELING.csv"
            write_csv(global_standardized, global_qc_output)
            output_files.append(global_qc_output)

        return {
            "split_pairs": split_pairs,
            "standardized_generated": False,
            "qc_output": qc_output,
            "output_files": output_files,
            "n_features_before_filtering": len(feature_cols),
            "n_features_after_missing_constant_qc": len([col for col in qc_df.columns if col not in metadata_cols]),
            "retained_feature_counts": retained_feature_counts,
        }

    standardized_dir.mkdir(parents=True, exist_ok=True)
    for seed, train_path, test_path in split_pairs:
        train_ids = _load_split_ids(train_path)
        test_ids = _load_split_ids(test_path)
        train_df = raw_df.loc[raw_df["molecule_id"].astype(str).isin(train_ids.astype(str))].copy()
        test_df = raw_df.loc[raw_df["molecule_id"].astype(str).isin(test_ids.astype(str))].copy()

        train_features, test_features, filtering_log, retained_features = _apply_training_filters(
            train_df[feature_cols],
            test_df[feature_cols],
            threshold=correlation_threshold,
        )

        train_output = standardized_dir / f"rdkit_descriptors_standardized_seed{seed}_train.csv"
        test_output = standardized_dir / f"rdkit_descriptors_standardized_seed{seed}_test.csv"
        scaler_output = processed_dir / f"feature_scaler_parameters_seed{seed}.csv"
        retained_output = processed_dir / f"rdkit_retained_features_seed{seed}.csv"
        filter_log_output = reports_dir / f"rdkit_feature_filtering_log_seed{seed}.csv"

        if retained_features:
            scaler = StandardScaler()
            train_scaled = scaler.fit_transform(train_features)
            test_scaled = scaler.transform(test_features)
            train_export = pd.concat(
                [
                    train_df[metadata_cols].reset_index(drop=True),
                    pd.DataFrame(train_scaled, columns=retained_features),
                ],
                axis=1,
            )
            test_export = pd.concat(
                [
                    test_df[metadata_cols].reset_index(drop=True),
                    pd.DataFrame(test_scaled, columns=retained_features),
                ],
                axis=1,
            )
            scaler_df = pd.DataFrame(
                {
                    "feature": retained_features,
                    "mean": scaler.mean_,
                    "scale": scaler.scale_,
                    "variance": scaler.var_,
                    "n_train": len(train_df),
                }
            )
        else:
            train_export = train_df[metadata_cols].reset_index(drop=True).copy()
            test_export = test_df[metadata_cols].reset_index(drop=True).copy()
            scaler_df = pd.DataFrame(columns=["feature", "mean", "scale", "variance", "n_train"])
        retained_df = pd.DataFrame({"feature": retained_features})

        write_csv(train_export, train_output)
        write_csv(test_export, test_output)
        write_csv(scaler_df, scaler_output)
        write_csv(retained_df, retained_output)
        write_csv(filtering_log, filter_log_output)

        output_files.extend([train_output, test_output, scaler_output, retained_output, filter_log_output])
        retained_feature_counts[seed] = len(retained_features)
        generated_standardized = True
        LOGGER.info(
            "Standardized RDKit descriptors for seed %s with %s retained features.",
            seed,
            len(retained_features),
        )

    return {
        "split_pairs": split_pairs,
        "standardized_generated": generated_standardized,
        "qc_output": qc_output,
        "output_files": output_files,
        "n_features_before_filtering": len(feature_cols),
        "n_features_after_missing_constant_qc": len([col for col in qc_df.columns if col not in metadata_cols]),
        "retained_feature_counts": retained_feature_counts,
    }
