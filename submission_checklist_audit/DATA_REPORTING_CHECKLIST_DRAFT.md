## Checklist item: Study design and reporting scope
Status: READY
Draft response:
Repository text explicitly describes benchmark prediction, structure-disjoint validation, external sensory-panel validation, and receptor-consistency analysis, and states the work is not a direct receptor mechanism proof. This appears supportable from the current repository snapshot.
Evidence in repository:
- `README.md`
- `docs/workflow_overview.md`
Remaining risk:
No major repository-level risk identified from the current audit.

## Checklist item: Dataset sources and access restrictions
Status: READY
Draft response:
Source inventory and third-party restriction note exist; PubChem and M2OR are mentioned in docs; label harmonization evidence exists in data/interim/label_harmonization_table.csv. This appears supportable from the current repository snapshot.
Evidence in repository:
- `data/raw_manifest/source_inventory.csv`
- `data/raw_manifest/third_party_source_restrictions.md`
- `docs/data_availability.md`
Remaining risk:
CHECK THIRD-PARTY RESTRICTION

## Checklist item: Sample size and dataset dimensions
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Current repository repeatedly states 4,495 benchmark molecules, 112 labels, 87 external molecules, 60 pairs, 43 panel labels, and 13 trained assessors. This item remains subject to author confirmation.
Evidence in repository:
- `README.md`
- `docs/workflow_overview.md`
- `docs/sensory_panel_metadata.md`
- `data/receptor_annotations/m2or_matched_molecules.csv`
Remaining risk:
Requested target counts 4,403 and 90 do not match observed release files.

## Checklist item: Inclusion/exclusion criteria
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Conflict-removal and M2OR exclusion evidence exists; workflow text mentions unresolved identifier-structure conflicts and QM exclusion logic, but qm_calculations logs are absent. This item remains subject to author confirmation.
Evidence in repository:
- `data/interim/removed_or_conflicting_entries.csv`
- `data/receptor_annotations/receptor_matching_exclusion_log.csv`
- `docs/workflow_overview.md`
Remaining risk:
Missing qm_calculations failure and imaginary-frequency log files.

## Checklist item: Data preprocessing and label harmonization
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
CAS/SMILES harmonization and SHAP ranking evidence exists. ECFP4 and RDKit assets exist. QM, RDKit+QM, scaler-parameter, correlation-filter, and ADCH public files requested by checklist are not present at the expected paths. This item remains subject to author confirmation.
Evidence in repository:
- `data/interim/cas_smiles_mapping.csv`
- `data/interim/molecule_records_canonicalized.csv`
- `data/interim/label_harmonization_table.csv`
- `results/shap/shap_feature_ranking.csv`
Remaining risk:
Partial release only.

## Checklist item: QM descriptor calculation and exclusion criteria
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Conflict-removal and M2OR exclusion evidence exists; workflow text mentions unresolved identifier-structure conflicts and QM exclusion logic, but qm_calculations logs are absent. This item remains subject to author confirmation.
Evidence in repository:
- `data/interim/removed_or_conflicting_entries.csv`
- `data/receptor_annotations/receptor_matching_exclusion_log.csv`
- `docs/workflow_overview.md`
Remaining risk:
Missing qm_calculations failure and imaginary-frequency log files.

## Checklist item: Molecular descriptor construction
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
CAS/SMILES harmonization and SHAP ranking evidence exists. ECFP4 and RDKit assets exist. QM, RDKit+QM, scaler-parameter, correlation-filter, and ADCH public files requested by checklist are not present at the expected paths. This item remains subject to author confirmation.
Evidence in repository:
- `data/interim/cas_smiles_mapping.csv`
- `data/interim/molecule_records_canonicalized.csv`
- `data/interim/label_harmonization_table.csv`
- `results/shap/shap_feature_ranking.csv`
Remaining risk:
Partial release only.

## Checklist item: Train/test split and leakage prevention
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Structure-disjoint train/test files exist. Leakage check result: cas overlap=0; smiles overlap=0. Test file exists. Workflow doc says scaling uses training-set parameters only. This item remains subject to author confirmation.
Evidence in repository:
- `data/splits/structure_disjoint_splits/structure_disjoint_train.csv`
- `data/splits/structure_disjoint_splits/structure_disjoint_test.csv`
- `tests/test_split_no_leakage.py`
- `docs/workflow_overview.md`
Remaining risk:
Iterative stratified split directory missing; no released feature_scaler_parameters.csv.

## Checklist item: Model training and hyperparameter reporting
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Configs list model families and random seeds 42-46; graph config covers GCN/ADCH variants. Multiple hyperparameter fields remain TODO, and trained weights/model cards are not present. This item remains subject to author confirmation.
Evidence in repository:
- `configs/model_hyperparameters.yaml`
- `configs/graph_model_settings.yaml`
- `scripts/`
- `results/benchmark_metrics/model_summary_json/`
Remaining risk:
Hyperparameters incomplete in released config; no weight files found.

## Checklist item: Evaluation metrics
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Benchmark AUC tables and receptor-consistency statistics exist. Graph settings mention Lin's concordance coefficient, Jaccard, and kNN macro-AUC, but the expected external_validation result files are not present. This item remains subject to author confirmation.
Evidence in repository:
- `results/benchmark_metrics/micro_auc_by_model_representation.csv`
- `results/benchmark_metrics/macro_auc_by_model_representation.csv`
- `results/benchmark_metrics/labelwise_delta_auc.csv`
- `results/receptor_consistency/receptor_consistency_statistics.csv`
- `configs/graph_model_settings.yaml`
Remaining risk:
No direct uncertainty / confidence interval artifact identified.

## Checklist item: External sensory panel
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Metadata explicitly states 13 trained assessors, 43-label vocabulary, knee-point Top-K extraction, aggregated anonymized release, and pairwise Jaccard summary. This item remains subject to author confirmation.
Evidence in repository:
- `docs/sensory_panel_metadata.md`
- `data/external_panel/aggregated_panel_ratings.csv`
- `data/external_panel/panel_derived_topk_labels.csv`
- `data/external_panel/pairwise_panel_jaccard_similarity.csv`
Remaining risk:
No explicit ethics approval, waiver, or informed consent statement found.

## Checklist item: Sensory panel ethics / consent / anonymization
Status: MISSING
Draft response:
Sensory panel involvement is documented, but no explicit ethics approval, waiver, informed consent, competing-interest, or funding statement was located in the scanned repository files. This item remains subject to author confirmation.
Evidence in repository:
- `docs/sensory_panel_metadata.md`
- `README.md`
Remaining risk:
HIGH PRIORITY MISSING ITEM

## Checklist item: Public receptor annotation matching
Status: READY
Draft response:
Processed M2OR matching outputs and exclusion log exist, and workflow text frames the analysis as orthogonal consistency rather than direct mechanistic proof. This appears supportable from the current repository snapshot.
Evidence in repository:
- `data/receptor_annotations/m2or_matched_molecules.csv`
- `data/receptor_annotations/receptor_response_vectors.csv`
- `data/receptor_annotations/receptor_pairwise_concordance.csv`
- `data/receptor_annotations/receptor_concordant_divergent_pairs.csv`
- `docs/receptor_matching_rules.md`
- `scripts/07_receptor_consistency_analysis.py`
Remaining risk:
No major repository-level risk identified from the current audit.

## Checklist item: Figure source data
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Fig. 1-4 source-data workbooks exist. No figures/ directory is present, supplementary_tables_source_data.xlsx is missing, and make_fig1.py to make_fig4.py were not found. This item remains subject to author confirmation.
Evidence in repository:
- `data/source_data/fig1_source_data.xlsx`
- `data/source_data/fig2_source_data.xlsx`
- `data/source_data/fig3_source_data.xlsx`
- `data/source_data/fig4_source_data.xlsx`
Remaining risk:
NEEDS MANUAL CHECK

## Checklist item: Code availability
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Core code and metadata files exist. There is no clear 00-09 one-click pipeline series, and version pinning is limited. This item remains subject to author confirmation.
Evidence in repository:
- `README.md`
- `requirements.txt`
- `environment.yml`
- `pyproject.toml`
- `scripts/`
- `src/qmodor/`
- `notebooks/`
- `tests/`
- `CITATION.cff`
- `.zenodo.json`
Remaining risk:
One-click reproduction script not found.

## Checklist item: Data availability
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Availability docs and processed data directories exist. .zenodo.json exists but related_identifiers is empty and no DOI is present. This item remains subject to author confirmation.
Evidence in repository:
- `docs/data_availability.md`
- `docs/data_dictionary.md`
- `data/raw_manifest/third_party_source_restrictions.md`
- `data/processed/`
- `data/source_data/`
Remaining risk:
Zenodo DOI or release tag missing.

## Checklist item: Software versions and computational environment
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Open-source package files exist; workflow doc names Gaussian 16W, Multiwfn 3.8 dev, and GaussView 6.0.16. The public environment files do not record all commercial-tool versions in machine-readable form, and CUDA/GPU usage is not clearly documented. This item remains subject to author confirmation.
Evidence in repository:
- `environment.yml`
- `requirements.txt`
- `docs/workflow_overview.md`
- `data/raw_manifest/third_party_source_restrictions.md`
Remaining risk:
Version granularity incomplete.

## Checklist item: Random seeds and reproducibility
Status: NEEDS AUTHOR CONFIRMATION
Draft response:
Configs list model families and random seeds 42-46; graph config covers GCN/ADCH variants. Multiple hyperparameter fields remain TODO, and trained weights/model cards are not present. This item remains subject to author confirmation.
Evidence in repository:
- `configs/model_hyperparameters.yaml`
- `configs/graph_model_settings.yaml`
- `scripts/`
- `results/benchmark_metrics/model_summary_json/`
Remaining risk:
Hyperparameters incomplete in released config; no weight files found.

## Checklist item: Competing interests / funding / author contributions
Status: MISSING
Draft response:
Sensory panel involvement is documented, but no explicit ethics approval, waiver, informed consent, competing-interest, or funding statement was located in the scanned repository files. This item remains subject to author confirmation.
Evidence in repository:
- `docs/sensory_panel_metadata.md`
- `README.md`
Remaining risk:
HIGH PRIORITY MISSING ITEM

## Checklist item: Third-party data redistribution restrictions
Status: READY
Draft response:
Source inventory and third-party restriction note exist; PubChem and M2OR are mentioned in docs; label harmonization evidence exists in data/interim/label_harmonization_table.csv. This appears supportable from the current repository snapshot.
Evidence in repository:
- `data/raw_manifest/source_inventory.csv`
- `data/raw_manifest/third_party_source_restrictions.md`
- `docs/data_availability.md`
Remaining risk:
CHECK THIRD-PARTY RESTRICTION
