"""Create traceable, non-destructive archival tables for the triplet experiment.

Inputs remain in ``三联体问题``.  This script only creates files in
``data/sensory_triplets`` and ``reports``.  It intentionally does not impute
responses, infer missing CAS/SMILES, or invent working concentrations.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import openpyxl
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "三联体问题"
OUT = ROOT / "data" / "sensory_triplets"
REPORTS = ROOT / "reports"

TRIPLET_CSV = SOURCE / "30组香原料三联体.csv"
TRIPLET_XLSX = SOURCE / "30组香原料三联体.xlsx"
RESPONSE_CSV = SOURCE / "Similarity_Sensory_Evaluation.csv"
RESPONSE_XLSX = SOURCE / "Similarity_Sensory_Evaluation.xlsx"
INTENSITY_CSV = SOURCE / "强度平衡.csv"
INTENSITY_XLSX = SOURCE / "强度平衡.xlsx"

CODE_MAP = {
    1: "A-P similar; A-S dissimilar",
    2: "A-S similar; A-P dissimilar",
    3: "Both similar",
    4: "Both dissimilar",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_triplets_csv(path: Path) -> pd.DataFrame:
    """Read the two-header-row source without treating group labels as fields."""
    rows = list(csv.reader(path.open("r", encoding="utf-8-sig", newline="")))
    required = {"序号", "编号", "CAS", "SMILES", "编号2", "CAS2", "SMILES2", "编号3", "CAS3", "SMILES3"}
    header_index = next(i for i, row in enumerate(rows[:10]) if required.issubset(set(row)))
    return pd.DataFrame(rows[header_index + 1 :], columns=rows[header_index]).replace({"": pd.NA})


def worksheet_values(path: Path) -> list[list[object]]:
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[wb.sheetnames[0]]
    return [[cell.value for cell in row] for row in ws.iter_rows()]


def normalized_rows(rows: list[list[object]]) -> list[list[str]]:
    def clean(value: object) -> str:
        if value is None:
            return ""
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value).strip()

    normalized = [[clean(value) for value in row] for row in rows]
    # Excel materializes the unused tail of the worksheet whereas CSV does not.
    # Trailing blank cells carry no table value and should not be a conflict.
    for row in normalized:
        while row and row[-1] == "":
            row.pop()
    return normalized


def source_equivalence(csv_path: Path, xlsx_path: Path, header_row: int = 0) -> tuple[bool, str]:
    csv_rows = list(csv.reader(csv_path.open("r", encoding="utf-8-sig", newline="")))
    xlsx_rows = worksheet_values(xlsx_path)
    a = normalized_rows(csv_rows[header_row:])
    b = normalized_rows(xlsx_rows[header_row:])
    if a == b:
        return True, f"CSV rows={len(a)}; XLSX rows={len(b)}; no value differences"
    for row_index, (csv_row, xlsx_row) in enumerate(zip(a, b), start=header_row + 1):
        if csv_row != xlsx_row:
            max_columns = max(len(csv_row), len(xlsx_row))
            for column_index in range(max_columns):
                left = csv_row[column_index] if column_index < len(csv_row) else "<absent>"
                right = xlsx_row[column_index] if column_index < len(xlsx_row) else "<absent>"
                if left != right:
                    return False, (f"CSV rows={len(a)}; XLSX rows={len(b)}; first difference "
                                   f"at source row {row_index}, column {column_index + 1}: "
                                   f"CSV={left!r}; XLSX={right!r}")
    return False, f"CSV rows={len(a)}; XLSX rows={len(b)}; row-count difference"


def source_differences(csv_path: Path, xlsx_path: Path, source_name: str, header_row: int = 0) -> list[dict[str, object]]:
    """Return every value-level CSV/XLSX difference, including row-length gaps."""
    csv_rows = normalized_rows(list(csv.reader(csv_path.open("r", encoding="utf-8-sig", newline="")))[header_row:])
    xlsx_rows = normalized_rows(worksheet_values(xlsx_path)[header_row:])
    differences = []
    for row_offset in range(max(len(csv_rows), len(xlsx_rows))):
        left_row = csv_rows[row_offset] if row_offset < len(csv_rows) else []
        right_row = xlsx_rows[row_offset] if row_offset < len(xlsx_rows) else []
        for col_offset in range(max(len(left_row), len(right_row))):
            left = left_row[col_offset] if col_offset < len(left_row) else "<absent>"
            right = right_row[col_offset] if col_offset < len(right_row) else "<absent>"
            if left != right:
                differences.append({
                    "Source_Table": source_name,
                    "Source_Row_1based": header_row + row_offset + 1,
                    "Source_Column_1based": col_offset + 1,
                    "CSV_Value": left,
                    "XLSX_Value": right,
                })
    return differences


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    required = [TRIPLET_CSV, TRIPLET_XLSX, RESPONSE_CSV, RESPONSE_XLSX, INTENSITY_CSV, INTENSITY_XLSX]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Required source file(s) absent: " + "; ".join(missing))

    # Definitions: retain only identity fields; all values are source-backed.
    raw_triplets = read_triplets_csv(TRIPLET_CSV)
    role_columns = {
        "A": ("编号", "名称", "CAS", "SMILES"),
        "S": ("编号2", "名称2", "CAS2", "SMILES2"),
        "P": ("编号3", "名称3", "CAS3", "SMILES3"),
    }
    definitions = pd.DataFrame({"Triplet_ID": raw_triplets["序号"].astype(int)})
    for role, (identifier, molecule, cas, smiles) in role_columns.items():
        definitions[f"{role}_Molecule_ID"] = raw_triplets[identifier]
        definitions[f"{role}_molecule"] = raw_triplets[molecule]
        definitions[f"{role}_CAS"] = raw_triplets[cas]
        definitions[f"{role}_SMILES"] = raw_triplets[smiles]
    definitions.to_csv(OUT / "triplet_definitions_30.csv", index=False, encoding="utf-8")

    # Responses: source uses evaluator1..13; release only anonymous Assessor_01..13.
    response_wide = pd.read_csv(RESPONSE_CSV, encoding="utf-8-sig")
    source_assessors = response_wide.columns[1:].tolist()
    assessor_map = {source: f"Assessor_{i:02d}" for i, source in enumerate(source_assessors, 1)}
    response_rows: list[dict[str, object]] = []
    for _, record in response_wide.iterrows():
        source_triplet = str(record.iloc[0])
        if not source_triplet.lower().startswith("combo"):
            raise ValueError(f"Unexpected triplet identifier: {source_triplet!r}")
        triplet_id = int(source_triplet[5:])
        for source_assessor, assessor_id in assessor_map.items():
            value = record[source_assessor]
            response_rows.append({
                "Triplet_ID": triplet_id,
                "Assessor_ID": assessor_id,
                "Response_Code": value,
                "Response_Category": CODE_MAP.get(int(value), pd.NA) if pd.notna(value) else pd.NA,
                "Source_Triplet_Label": source_triplet,
                "Source_Assessor_Column": source_assessor,
            })
    responses = pd.DataFrame(response_rows)
    responses.to_csv(OUT / "triplet_assessor_responses_390.csv", index=False, encoding="utf-8")
    pd.DataFrame(
        [{"Source_Assessor_Column": source, "Assessor_ID": anonymous} for source, anonymous in assessor_map.items()]
    ).to_csv(OUT / "assessor_anonymization_map.csv", index=False, encoding="utf-8")

    # Summary is derived only from the unmodified codes in the preceding table.
    valid = responses[responses["Response_Code"].isin(CODE_MAP)].copy()
    summary_rows = []
    for triplet_id, group in valid.groupby("Triplet_ID", sort=True):
        counts = Counter(group["Response_Code"].astype(int))
        n = len(group)
        ap = (counts[1] + counts[3]) / n if n else pd.NA
        a_s = (counts[2] + counts[3]) / n if n else pd.NA
        summary_rows.append({
            "Triplet_ID": triplet_id,
            "N_Assessors": n,
            "Count_AP_Similar": counts[1],
            "Count_AS_Similar": counts[2],
            "Count_Both_Similar": counts[3],
            "Count_Both_Dissimilar": counts[4],
            "Proportion_AP_Similar": ap,
            "Proportion_AS_Similar": a_s,
            "Perceptual_Similarity_Difference": ap - a_s if n else pd.NA,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "triplet_response_summary.csv", index=False, encoding="utf-8")

    # Intensity source contains one reported score per molecule/concentration,
    # not assessor-level data, averages, or selected working concentrations.
    intensity_wide = pd.read_csv(INTENSITY_CSV, encoding="utf-8-sig")
    molecule_lookup = {}
    for _, row in definitions.iterrows():
        for role in ("A", "S", "P"):
            molecule_lookup[row[f"{role}_Molecule_ID"]] = {
                "Molecule_Name": row[f"{role}_molecule"], "CAS": row[f"{role}_CAS"]
            }
    concentration_cols = [column for column in intensity_wide.columns if column.endswith("强度")]
    intensity_rows = []
    for _, row in intensity_wide.iterrows():
        molecule_id = row["分子编号"]
        info = molecule_lookup.get(molecule_id, {})
        for column in concentration_cols:
            intensity_rows.append({
                "Molecule_ID": molecule_id,
                "Molecule_Name_CN": row["中文名称"],
                "Molecule_Name_EN": row["标准英文名称"],
                "CAS": info.get("CAS", pd.NA),
                "Concentration": column.removesuffix("强度"),
                "Assessor_ID": pd.NA,
                "Intensity_Score": row[column],
                "Measure_Type": "reported_score_from_source_matrix",
                "Selected_Working_Concentration": pd.NA,
                "Source_Column": column,
            })
    intensity = pd.DataFrame(intensity_rows)
    intensity.to_csv(OUT / "triplet_intensity_balancing.csv", index=False, encoding="utf-8")

    # Checks: do not repair any anomaly; report all observed status fields.
    definition_ids = set(definitions["Triplet_ID"])
    response_ids = set(responses["Triplet_ID"])
    numeric_codes = pd.to_numeric(responses["Response_Code"], errors="coerce")
    valid_code_mask = numeric_codes.isin(CODE_MAP)
    duplicate_response_rows = int(responses.duplicated(["Triplet_ID", "Assessor_ID"]).sum())
    # The CSV's first grouping header uses literal "Unnamed:" labels where the
    # workbook has merged/blank cells. Compare the actual field header onward.
    triplet_csv_xlsx_equal, triplet_compare_detail = source_equivalence(TRIPLET_CSV, TRIPLET_XLSX, header_row=1)
    response_csv_xlsx_equal, response_compare_detail = source_equivalence(RESPONSE_CSV, RESPONSE_XLSX)
    intensity_csv_xlsx_equal, intensity_compare_detail = source_equivalence(INTENSITY_CSV, INTENSITY_XLSX)
    conflicts = (
        source_differences(TRIPLET_CSV, TRIPLET_XLSX, "triplet_definitions", header_row=1)
        + source_differences(RESPONSE_CSV, RESPONSE_XLSX, "sensory_responses")
        + source_differences(INTENSITY_CSV, INTENSITY_XLSX, "intensity_balancing")
    )
    pd.DataFrame(conflicts, columns=["Source_Table", "Source_Row_1based", "Source_Column_1based", "CSV_Value", "XLSX_Value"]).to_csv(
        REPORTS / "triplet_source_csv_xlsx_conflicts.csv", index=False, encoding="utf-8"
    )
    per_triplet = valid.groupby("Triplet_ID").size()
    source_molecules = set(pd.concat([definitions[f"{role}_Molecule_ID"] for role in ("A", "S", "P")]))
    intensity_molecules = set(intensity_wide["分子编号"])

    checks = {
        "triplet_count_is_30": len(definitions) == 30,
        "response_count_is_390": len(responses) == 390,
        "all_triplets_have_13_valid_responses": len(per_triplet) == 30 and bool((per_triplet == 13).all()),
        "all_response_codes_in_1_to_4": bool(valid_code_mask.all()),
        "response_triplet_ids_match_definition_ids": response_ids == definition_ids,
        "no_duplicate_triplet_assessor_records": duplicate_response_rows == 0,
        "each_summary_row_reconciles_to_n_assessors": bool((summary[["Count_AP_Similar", "Count_AS_Similar", "Count_Both_Similar", "Count_Both_Dissimilar"]].sum(axis=1) == summary["N_Assessors"]).all()),
        "triplet_csv_equals_xlsx": triplet_csv_xlsx_equal,
        "response_csv_equals_xlsx": response_csv_xlsx_equal,
        "intensity_csv_equals_xlsx": intensity_csv_xlsx_equal,
        "intensity_molecules_match_triplet_molecule_set": intensity_molecules == source_molecules,
        "intensity_long_row_count_is_255": len(intensity) == len(intensity_wide) * len(concentration_cols),
        "all_outputs_utf8_decodable": True,
    }
    for output in OUT.glob("*.csv"):
        output.read_text(encoding="utf-8")

    # The README documents both the evidence and limitations rather than
    # upgrading the source matrix into data that it does not contain.
    (OUT / "README.md").write_text(f"""# Structure-odor discontinuity triplet experiment

## Scope and provenance

This directory is a machine-readable, non-destructive archive of the triplet sensory experiment. It was generated from source files retained in `三联体问题/`:

- `30组香原料三联体.csv` / `.xlsx`: 30 A/S/P triplet definitions.
- `Similarity_Sensory_Evaluation.csv` / `.xlsx`: 30 triplets by 13 evaluators.
- `强度平衡.csv` / `.xlsx`: 51 molecules by five concentration columns.

The files were compared as CSV versus first worksheet of XLSX. Results: triplets={triplet_csv_xlsx_equal}; responses={response_csv_xlsx_equal}; intensity={intensity_csv_xlsx_equal}. The normalized archive uses CSV because `plot.py` directly reads that response file. This is a provenance choice, not an adjudication that CSV is scientifically authoritative; see `reports/triplet_source_csv_xlsx_conflicts.csv` for every unresolved cell difference. SHA-256 values and detailed checks are in `reports/triplet_data_audit_report.md`.

## Experiment design and response coding

Each triplet has an anchor (A), structural neighbour (S), and perceptual neighbour (P). The source sensory matrix contains 30 triplets × 13 evaluator columns = 390 cells. Source evaluator labels are replaced by `Assessor_01` through `Assessor_13`; the non-identifying source-column-to-release-ID mapping is retained in `assessor_anonymization_map.csv`.

The mapping is confirmed by the comments and legend in source `plot.py`, and corroborated by the text extracted from `Figure3B_sensory_response_matrix.pdf`: code 1 = A-P similar / A-S dissimilar; code 2 = A-S similar / A-P dissimilar; code 3 = both similar; code 4 = both dissimilar. The source spreadsheets do not themselves contain a prose code legend.

## Files

- `triplet_definitions_30.csv`: One source triplet per row. `Triplet_ID` is source sequence number; identity, CAS, and SMILES columns are retained separately for A, S, and P. No absent identifier has been inferred.
- `triplet_assessor_responses_390.csv`: One source matrix cell per row. `Response_Code` is unaltered; `Response_Category` is the verified mapping above. `Source_Triplet_Label` and `Source_Assessor_Column` preserve source traceability.
- `triplet_response_summary.csv`: Counts and derived proportions grouped from the response file. `Proportion_AP_Similar=(code1+code3)/N`; `Proportion_AS_Similar=(code2+code3)/N`; the difference is AP minus AS. This is a new transparent summary, not a replacement for a published figure statistic.
- `triplet_intensity_balancing.csv`: Long-form transcription of the intensity matrix. Every `Intensity_Score` is exactly one supplied cell. No `Assessor_ID`, replicate count, mean/SD label, or selected working concentration exists in the supplied source, so those fields remain empty and `Measure_Type` identifies the source-level granularity. Concentration is a source label (`0.5%` through `10.0%`), not a converted physical unit.

## Figure relationship

`plot.py` reads the source response CSV and creates `Figure3B_sensory_response_matrix.png` and `.pdf`; these are visualization products, not a second raw data source. Their historical `Figure3B` name is retained. In the current paper context, the triplet result is described as Figure 4; the filename alone does not establish the current figure number. The intensity table is the available input relevant to Supplementary Figure S2. No source mapping to Supplementary Figure S3 was found in this folder.

## Reproduction

Run `python scripts/audit_triplet_sensory_data.py` from repository root using an environment with pandas and openpyxl. The script neither edits source files nor trains models. See the audit report for checks and unresolved items.
""", encoding="utf-8")

    file_rows = []
    for path in sorted(SOURCE.iterdir()):
        if path.is_file():
            file_rows.append((path.name, path.stat().st_size, sha256(path)))
    status_lines = "\n".join(f"| {name} | {size} | `{digest}` |" for name, size, digest in file_rows)
    check_lines = "\n".join(f"| {name} | {'PASS' if passed else 'FAIL'} |" for name, passed in checks.items())
    failed = [name for name, passed in checks.items() if not passed]
    (REPORTS / "triplet_data_audit_report.md").write_text(f"""# Triplet sensory-data audit report

## Source inventory

All source files were found in `{SOURCE.relative_to(ROOT)}` and were read without modification. The project root has a `.git` directory, but `git status` did not recognize this directory as a valid Git work tree during this audit; no Git operation was performed.

| Source file | Bytes | SHA-256 |
| --- | ---: | --- |
{status_lines}

## Source content and relationship

| Source | Actual content | Role |
| --- | --- | --- |
| `30组香原料三联体.csv/.xlsx` | 2 header rows + {len(raw_triplets)} triplet rows; 29 columns | Raw triplet definition/selection table. First header is a grouping row. |
| `Similarity_Sensory_Evaluation.csv/.xlsx` | {len(response_wide)} triplets × {len(source_assessors)} evaluator columns | Raw sensory response matrix. |
| `强度平衡.csv/.xlsx` | {len(intensity_wide)} molecule rows × {len(concentration_cols)} concentration-score columns | Reported intensity score matrix; no evaluator/replicate/selection columns. |
| `plot.py` | Reads sensory CSV and generates response matrix PNG/PDF | Plotting product; also documents code legend. |
| `make_triplet_2d_excel.py` | Reads triplet CSV and can create structure-image Excel file | Optional visualization utility, not used to derive archive tables. It has default relative file names and needs RDKit/xlsxwriter. |
| `xlsx_to_csv_converter.py` | Converts every first-sheet XLSX in its current directory to CSV | Generic conversion utility; it overwrites matching CSV targets if run in source folder. Not used. |
| `Figure3B_sensory_response_matrix.png/.pdf` | Rendered 30×13 response matrix | Plot products generated by `plot.py`; PDF has one page and its legend corroborates code mapping. |

## Encoding verification

`plot.py` explicitly validates codes {{1,2,3,4}} and defines: 1 A-P similar/A-S dissimilar; 2 A-S similar/A-P dissimilar; 3 both similar; 4 both dissimilar. The generated PDF legend contains the same text (with extraction-character substitution for en dashes). This mapping was therefore used only to add `Response_Category`; numeric source codes remain unaltered.

## Automated checks

| Check | Status |
| --- | --- |
{check_lines}

CSV/XLSX comparison details: triplets ({triplet_compare_detail}); responses ({response_compare_detail}); intensity ({intensity_compare_detail}).

`reports/triplet_source_csv_xlsx_conflicts.csv` contains {len(conflicts)} cell-level differences. The archive therefore uses the CSV versions, because `plot.py` reads `Similarity_Sensory_Evaluation.csv`; this does not resolve the contradiction. The response conflict is scientifically material and must be adjudicated before treating CSV/XLSX as interchangeable source files.

Observed counts: {len(definitions)} triplets, {len(responses)} response cells, {len(valid)} valid response codes, {len(source_molecules)} distinct source molecules across A/S/P, and {len(intensity_wide)} intensity-matrix molecules. Response-code counts: {dict(sorted(Counter(valid['Response_Code'].astype(int)).items()))}. Missing response cells: {int(responses['Response_Code'].isna().sum())}; duplicate triplet-assessor keys: {duplicate_response_rows}; illegal/non-numeric codes: {int((~valid_code_mask).sum())}.

## Output mapping

| Output | Derived from | Content |
| --- | --- | --- |
| `data/sensory_triplets/triplet_definitions_30.csv` | triplet source CSV | Identity/CAS/SMILES fields for all A/S/P rows. |
| `data/sensory_triplets/triplet_assessor_responses_390.csv` | response source CSV | Lossless long-form expansion of the 30×13 response matrix plus non-identifying trace columns. |
| `data/sensory_triplets/triplet_response_summary.csv` | normalized response table | Per-triplet response counts and documented derived proportions. |
| `data/sensory_triplets/triplet_intensity_balancing.csv` | intensity source CSV + definition lookup | Long-form score matrix with CAS only when traceable through triplet definitions. |
| `data/sensory_triplets/README.md` | sources and this audit | Release documentation and figure/history notes. |
| `reports/triplet_source_csv_xlsx_conflicts.csv` | direct CSV/XLSX comparison | Every unequal source cell; no value was changed to resolve it. |

## Unresolved items requiring human confirmation

- The supplied intensity table has no assessor-level observations, replicate counts, aggregation method, or selected working concentration. It cannot independently reproduce a raw-level Supplementary Figure S2 analysis or establish final concentrations.
- The source folder does not state a mapping to current Figure 4 or Supplementary Figure S3. `Figure3B` is treated as a historical filename only.
- No direct experimental protocol or source worksheet prose independently describes the 1-4 legend; verification rests on the plotting script and its generated PDF.
- The project repository Git metadata appears incomplete or nonstandard from this shell (`git status` failed), so a baseline dirty-worktree status could not be established.
- No names were found in the sensory response source. The released assessor IDs are anonymized column aliases.

## Concise file list

The candidate public tables are `triplet_definitions_30.csv`, `triplet_assessor_responses_390.csv`, `triplet_response_summary.csv`, `triplet_intensity_balancing.csv`, and the accompanying `README.md`. Review the unresolved items before public release.
""", encoding="utf-8")

    (OUT / "audit_metadata.json").write_text(json.dumps({
        "source_file_sha256": {name: digest for name, _, digest in file_rows},
        "checks": checks,
        "response_code_counts": {str(k): v for k, v in sorted(Counter(valid["Response_Code"].astype(int)).items())},
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"checks": checks, "failed_checks": failed}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
