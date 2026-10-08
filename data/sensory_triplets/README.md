# Structure-odor discontinuity triplet experiment

## Scope and provenance

This directory is a machine-readable, non-destructive archive of the triplet sensory experiment. It was generated from source files retained in `三联体问题/`:

- `30组香原料三联体.csv` / `.xlsx`: 30 A/S/P triplet definitions.
- `Similarity_Sensory_Evaluation.csv` / `.xlsx`: 30 triplets by 13 evaluators.
- `强度平衡.csv` / `.xlsx`: 51 molecules by five concentration columns.

The files were compared as CSV versus first worksheet of XLSX. Results: triplets=False; responses=False; intensity=True. The normalized archive uses CSV because `plot.py` directly reads that response file. This is a provenance choice, not an adjudication that CSV is scientifically authoritative; see `reports/triplet_source_csv_xlsx_conflicts.csv` for every unresolved cell difference. SHA-256 values and detailed checks are in `reports/triplet_data_audit_report.md`.

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
