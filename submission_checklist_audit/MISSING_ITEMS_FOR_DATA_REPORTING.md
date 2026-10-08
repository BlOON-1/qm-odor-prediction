HIGH PRIORITY:
- 问题: Ethics / informed consent statement not found
- 当前证据: docs/sensory_panel_metadata.md documents anonymized aggregated release only; no explicit ethics approval / waiver / consent statement located.
- 为什么是风险: Human sensory assessor work often triggers editorial compliance questions.
- 建议作者补充什么: Provide the exact ethics approval, waiver, or informed-consent wording used for the panel study.
- 不要自行补写的内容: Do not invent approval numbers, committee names, or consent wording.
- 问题: Zenodo DOI or release tag missing
- 当前证据: .zenodo.json exists but has no DOI in related_identifiers and no release tag evidence was found in the repository snapshot.
- 为什么是风险: Editors may require a persistent archive identifier for code/data availability.
- 建议作者补充什么: Mint or confirm the archival DOI / release tag and align the checklist wording with that identifier.
- 不要自行补写的内容: Do not guess a DOI or release URL.
- 问题: Leakage prevention evidence incomplete
- 当前证据: Structure-disjoint files and test file exist, but iterative_stratified_splits files and feature_scaler_parameters.csv are missing; training-only scaling is stated in docs but not fully auditable from files alone.
- 为什么是风险: Split provenance and leakage controls are common editorial review points for ML studies.
- 建议作者补充什么: Provide the missing iterative split artifacts or a clear author statement describing where they are archived and how scaling was fit only on training data.
- 不要自行补写的内容: Do not claim train-only scaling was verified if the file evidence is absent.
- 问题: Third-party restriction statement must be carried into submission text
- 当前证据: third_party_source_restrictions.md exists and clearly limits redistribution of TGSC, Leffingwell, M2OR, PubChem-linked materials, and QM intermediates.
- 为什么是风险: Checklist answers that overstate public redistributability could create compliance problems.
- 建议作者补充什么: Include a short restriction statement in the Data Reporting Checklist and data availability responses.
- 不要自行补写的内容: Do not paste raw third-party source tables into the submission materials.

MEDIUM PRIORITY:
- 问题: Hyperparameter detail remains incomplete
- 当前证据: configs/model_hyperparameters.yaml contains many TODO placeholders for penalty, C, gamma, n_estimators, learning_rate, epochs, batch_size, and related details.
- 为什么是风险: Reviewers may ask how final model settings were chosen and whether test data influenced tuning.
- 建议作者补充什么: Add the exact final hyperparameters or prepare a checklist note stating where they are archived.
- 不要自行补写的内容: Do not fill TODO fields with guessed values.
- 问题: Software version detail is incomplete
- 当前证据: Workflow docs name Gaussian 16W, Multiwfn 3.8 dev, and GaussView 6.0.16, but environment files do not fully encode all versions or GPU/CUDA details.
- 为什么是风险: Version ambiguity can weaken reproducibility claims.
- 建议作者补充什么: Provide final software version strings and whether GPU/CUDA was used for any reported training.
- 不要自行补写的内容: Do not assume versions from memory or local machine defaults.
- 问题: Supplementary tables source-data workbook missing
- 当前证据: data/source_data/fig1_source_data.xlsx through fig4_source_data.xlsx exist, but supplementary_tables_source_data.xlsx was not found.
- 为什么是风险: If the journal expects a checklist-level mapping for supplementary tables, this will look incomplete.
- 建议作者补充什么: Confirm whether supplementary table source data are included elsewhere or need a dedicated workbook.
- 不要自行补写的内容: Do not claim the workbook exists when it is absent.
- 问题: Figure traceability remains partly manual
- 当前证据: No figures/ directory and no make_fig1.py to make_fig4.py scripts were found in the scanned repository snapshot.
- 为什么是风险: Editors may ask how figure PDFs/source data connect to result tables.
- 建议作者补充什么: Prepare a short manual traceability note from each figure workbook to its source result file(s).
- 不要自行补写的内容: Do not assert full numeric traceability without a manual check.

LOW PRIORITY:
- 问题: Repository narrative could more explicitly explain release-vs-manuscript count differences
- 当前证据: Observed files say 4,495 benchmark molecules and 87 external molecules, while the checklist request references 4,403 and 90.
- 为什么是风险: Count mismatches can confuse editors even when scientifically explainable.
- 建议作者补充什么: Add a concise author confirmation note explaining which counts belong in the checklist and why.
- 不要自行补写的内容: Do not silently substitute one set of counts for another.
