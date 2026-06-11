"""Unified entry point for RDKit descriptor and ECFP4 generation."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT_DEFAULT = SCRIPT_PATH.parent.parent
SRC_DIR = PROJECT_ROOT_DEFAULT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qmodor.descriptors.common import DependencyError, load_benchmark_table, write_text
from qmodor.descriptors.compute_ecfp4 import compute_ecfp4
from qmodor.descriptors.compute_rdkit_descriptors import compute_rdkit_descriptors
from qmodor.descriptors.standardize_rdkit_descriptors import standardize_rdkit_descriptors


def parse_args() -> argparse.Namespace:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(description="Compute RDKit descriptors and ECFP4 fingerprints.")
    parser.add_argument("--project-root", default=str(PROJECT_ROOT_DEFAULT))
    parser.add_argument("--input", default=None)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n-bits", type=int, default=1024)
    parser.add_argument("--correlation-threshold", type=float, default=0.9)
    parser.add_argument("--allow-global-qc-standardization", action="store_true")
    return parser.parse_args()


def build_generation_report(
    project_root: Path,
    rdkit_result: dict,
    ecfp4_result: dict,
    standardization_result: dict,
) -> Path:
    """Write a markdown report summarizing RDKit and ECFP4 generation outputs."""

    raw_shape = rdkit_result["raw_df"].shape
    ecfp4_shape = ecfp4_result["matrix_shape"]
    before_filter = len(rdkit_result["descriptor_names_before_filtering"])
    after_filter = standardization_result["n_features_after_missing_constant_qc"]
    split_found = bool(standardization_result["split_pairs"])
    standardized_generated = bool(standardization_result["standardized_generated"])

    files_generated = sorted(
        {
            str(path.relative_to(project_root)).replace("\\", "/")
            for path in (
                rdkit_result["output_files"]
                + ecfp4_result["output_files"]
                + standardization_result["output_files"]
            )
        }
    )
    files_not_generated = []
    if not split_found:
        files_not_generated.extend(
            [
                "data/processed/standardized/rdkit_descriptors_standardized_seed42_train.csv and related seed-specific outputs: split files were not available.",
                "data/processed/feature_scaler_parameters_seed*.csv: scaler fitting on all molecules was not performed to avoid data leakage.",
            ]
        )

    warnings = []
    if len(rdkit_result["invalid_df"]) > 0:
        warnings.append(
            f"Invalid or missing SMILES were encountered for {len(rdkit_result['invalid_df'])} molecules."
        )
    if not split_found:
        warnings.append(
            "Train/test split files were not found, so finalized standardized RDKit descriptor matrices were not generated."
        )

    report_path = project_root / "reports" / "rdkit_ecfp4_generation_report.md"
    files_generated.append("reports/rdkit_ecfp4_generation_report.md")

    report_lines = [
        "# RDKit and ECFP4 Generation Report",
        "",
        "## 1. Input table used",
        f"- `{rdkit_result['input_info'].path.relative_to(project_root).as_posix()}`",
        "",
        "## 2. Molecule count",
        f"- total molecules: `{len(rdkit_result['input_info'].df)}`",
        "",
        "## 3. Valid molecule count",
        f"- valid molecules: `{len(rdkit_result['raw_df'])}`",
        "",
        "## 4. Invalid SMILES count",
        f"- invalid SMILES: `{len(rdkit_result['invalid_df'])}`",
        "",
        "## 5. RDKit descriptor count before filtering",
        f"- descriptor count before filtering: `{before_filter}`",
        "",
        "## 6. RDKit descriptor count after missing/constant filtering",
        f"- descriptor count after missing/constant filtering: `{after_filter}`",
        "",
        "## 7. Whether split files were found",
        f"- split files found: `{split_found}`",
        "",
        "## 8. Whether standardized descriptors were generated",
        f"- standardized descriptors generated: `{standardized_generated}`",
        "",
        "## 9. ECFP4 settings: radius=2, nBits=1024",
        f"- radius: `{ecfp4_result['metadata']['radius']}`",
        f"- nBits: `{ecfp4_result['metadata']['nBits']}`",
        f"- matrix shape: `{ecfp4_shape}`",
        "",
        "## 10. Files generated",
    ]
    report_lines.extend([f"- `{path}`" for path in files_generated])
    report_lines.extend(["", "## 11. Files not generated and why"])
    if files_not_generated:
        report_lines.extend([f"- {item}" for item in files_not_generated])
    else:
        report_lines.append("- None.")
    report_lines.extend(["", "## 12. Critical warnings"])
    if warnings:
        report_lines.extend([f"- {warning}" for warning in warnings])
    else:
        report_lines.append("- None.")

    write_text("\n".join(report_lines) + "\n", report_path)
    return report_path


def main() -> int:
    """Run RDKit descriptor generation, ECFP4 generation, and conditional standardization."""

    args = parse_args()
    project_root = Path(args.project_root).resolve()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    try:
        input_info = load_benchmark_table(project_root, args.input)
        rdkit_result = compute_rdkit_descriptors(project_root=project_root, input_path=args.input)
        ecfp4_result = compute_ecfp4(
            project_root=project_root,
            input_path=args.input,
            radius=args.radius,
            n_bits=args.n_bits,
        )
        standardization_result = standardize_rdkit_descriptors(
            project_root=project_root,
            correlation_threshold=args.correlation_threshold,
            allow_global_qc_standardization=args.allow_global_qc_standardization,
        )
        report_path = build_generation_report(project_root, rdkit_result, ecfp4_result, standardization_result)
    except (DependencyError, FileNotFoundError, ValueError) as exc:
        logging.error(str(exc))
        return 1

    output_files = sorted(
        {
            str(path.relative_to(project_root)).replace("\\", "/")
            for path in (
                rdkit_result["output_files"]
                + ecfp4_result["output_files"]
                + standardization_result["output_files"]
                + [report_path]
            )
        }
    )

    print(f"input file used: {input_info.path.relative_to(project_root).as_posix()}")
    print(f"total molecules: {len(input_info.df)}")
    print(f"valid RDKit molecules: {len(rdkit_result['raw_df'])}")
    print(f"invalid SMILES count: {len(rdkit_result['invalid_df'])}")
    print(f"RDKit raw descriptor shape: {rdkit_result['raw_df'].shape}")
    print(f"ECFP4 matrix shape: {ecfp4_result['matrix_shape']}")
    print(f"whether split files were found: {bool(standardization_result['split_pairs'])}")
    print(f"whether standardized descriptors were generated: {bool(standardization_result['standardized_generated'])}")
    print("output files generated:")
    for path in output_files:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
