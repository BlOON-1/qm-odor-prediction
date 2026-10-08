# Files to fix, exclude, or document

## Missing planned assets

| status | asset | evidence | recommendation | affects reproduction |
|---|---|---|---|---|
| MISSING (critical) | QM descriptor matrix; ADCH charge/features; RDKit+QM fusion matrix | No matching released files found | Export from verified source with IDs and QC/failure log | Yes |
| MISSING (critical) | Per-molecule prediction probabilities, trained weights/seeds, SHAP values | No matching released files found | Release verified outputs or limit claims | Yes |
| MISSING (critical) | External-panel assessor-level records | Only aggregated panel values found | Release anonymized raw ratings or state aggregate-only limitation | Yes |
| MISSING (important) | Docking validation source records for 8F76, 8UXY, 9WG4 | No matching validation inputs/results found | Add source/RMSD/pocket-recovery records | Yes |
| MISSING (optional) | Original triplet Excel source | Derived CSV exists, but source workbook absent | Add licensed/anonymized source workbook or provenance hash | No |

## Existing files not approved for direct upload

| status | file | evidence | recommendation | affects reproduction |
|---|---|---|---|---|
| UNVERIFIED | `.gitignore` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `.zenodo.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `CITATION.cff` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/dataset.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/descriptor_sets.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/figure_settings.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/graph_model_settings.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/model_hyperparameters.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/qm_workflow.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `configs/split_settings.yaml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/docking/final_results/all_failures.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/binding_mode_comparison.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/cluster_summary.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/comparison_discovery_status.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/comparison_eligibility.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/comparison_skips.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/cutoff_sensitivity.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/energy_summary.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/FINAL_REPORT.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/FINAL_REPORT.txt` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/input_manifest.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/interaction_summary.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/output_manifest.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/PIPELINE_STATUS.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/quality_control.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/receptor_ligand_groups.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/receptor_ligand_index.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/residue_interaction_frequency.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/RESULT_FILES_README.txt` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/final_results/software_versions.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/input_manifest_portable.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/receptor_ligand_index_portable.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/source_input_manifest.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/source_pipeline_jobs.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/source_software_versions.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/clustering.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/comparison.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/complexes.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/docking.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/extraction.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/plip.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/post_analysis.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/preflight.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/receptor_index.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/report.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/metadata/stage_status/validation.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/halogen_bonds.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/hydrogen_bonds.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/hydrophobic_contacts.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/interaction_summary.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/metal_interactions.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/pi_cation.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/pi_stacking.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/plip_execution_status.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/plip_parse_failures.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/residue_interaction_fingerprint.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/residue_interaction_frequency.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/docking/plip_details/salt_bridges.tsv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/external_panel/aggregated_panel_ratings.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/external_panel/aggregated_panel_ratings_wide_mean.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/external_panel/external_molecules_observed.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/external_panel/panel_derived_topk_labels.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/interim/cas_smiles_mapping.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/interim/label_harmonization_table.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/interim/molecule_records_canonicalized.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/interim/removed_or_conflicting_entries.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/processed/benchmark_molecules_observed.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/processed/odor_label_matrix_observed.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/processed/odor_label_vocabulary_observed.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/processed/rdkit_descriptors_raw.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/processed/rdkit_descriptors_unstandardized_qc.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/raw_manifest/source_inventory.csv` | Not validated sufficiently for PASS | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/raw_manifest/third_party_source_restrictions.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/receptor_annotations/m2or_matched_molecules.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/m2or_matched_molecules_all_species.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/m2or_pair_or_raw.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_concordant_divergent_pairs.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_concordant_divergent_pairs_all_species.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_matching_exclusion_log.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_pairwise_concordance.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_pairwise_concordance_all_species.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_response_long.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_response_long_all_species.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_response_vectors.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/receptor_annotations/receptor_response_vectors_all_species.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/sensory_triplets/assessor_anonymization_map.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/sensory_triplets/audit_metadata.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/sensory_triplets/README.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/sensory_triplets/triplet_intensity_balancing.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/source_data/fig1_source_data.xlsx` | Not validated sufficiently for PASS | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/source_data/fig2_source_data.xlsx` | Not validated sufficiently for PASS | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/source_data/fig3_source_data.xlsx` | Not validated sufficiently for PASS | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/source_data/fig4_source_data.xlsx` | Not validated sufficiently for PASS | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `data/splits/structure_disjoint_splits/fold_01_test.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/splits/structure_disjoint_splits/fold_01_train.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/splits/structure_disjoint_splits/structure_disjoint_split_audit.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/splits/structure_disjoint_splits/structure_disjoint_test.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `data/splits/structure_disjoint_splits/structure_disjoint_train.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `docs/data_availability.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/data_dictionary.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/qm_descriptor_definitions.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/receptor_matching_rules.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/reproducibility_checklist.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/sensory_panel_metadata.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `docs/workflow_overview.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `environment.yml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `LICENSE` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `notebooks/01_dataset_overview.ipynb` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `notebooks/02_descriptor_quality_control.ipynb` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `notebooks/03_benchmark_result_summary.ipynb` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `notebooks/04_external_panel_validation.ipynb` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `notebooks/05_receptor_consistency_check.ipynb` | binary/non-tabular asset not structurally parsed | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `pyproject.toml` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `README.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `requirements.txt` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `results/benchmark_metrics/labelwise_auc_by_model_representation.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/labelwise_delta_auc.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/macro_auc_by_model_representation.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/micro_auc_by_model_representation.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/CNN__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/CNN__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/CNN__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/CNN__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/LogReg__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/LogReg__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/LogReg__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/LogReg__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/MLP__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/MLP__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/MLP__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/MLP__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/RF__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/RF__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/RF__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/RF__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/summary_json_manifest.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/SVM__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/SVM__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/SVM__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/SVM__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/XGBoost__ECFP4__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/XGBoost__QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/XGBoost__RDKit__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/benchmark_metrics/model_summary_json/XGBoost__RDKit_QM__summary.json` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/receptor_consistency/receptor_consistency_statistics.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/shap/shap_feature_ranking.csv` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `results/structure_disjoint/README_pending_manual_confirmation.md` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | Yes |
| UNVERIFIED | `scripts/01_compute_rdkit_features.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/03_make_splits.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/07_receptor_consistency_analysis.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/audit_and_stage_docking.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/audit_triplet_sensory_data.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/build_complexes.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/cluster_poses.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/compare_binding_modes.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/extract_dlg_poses.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/generate_final_report.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/parse_plip_results.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/pipeline_common.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/pipeline_config.sh` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/run_autodock4_batch_pipeline_64core.sh` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/run_plip_batch.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/docking/validate_pipeline_outputs.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/legacy_split_ecfp4_by_folds.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/prepare_available_data.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/release_integrity_audit.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `scripts/sanitize_release_provenance.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/__init__.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/descriptors/__init__.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/descriptors/common.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/descriptors/compute_ecfp4.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/descriptors/compute_rdkit_descriptors.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/descriptors/standardize_rdkit_descriptors.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/receptor/__init__.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/receptor/compute_receptor_descriptor_distances.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/receptor/compute_receptor_pair_concordance.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/receptor/prepare_m2or_receptor_annotations.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/splits/__init__.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/splits/legacy_split_ecfp4_by_folds.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `src/qmodor/splits/organize_structure_disjoint_split.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `tests/test_descriptor_shapes.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `tests/test_receptor_matching.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |
| UNVERIFIED | `tests/test_split_no_leakage.py` | readable, but no file-specific provenance/cross-source validation | Review source/provenance; redact or exclude where applicable. | No |

Do not repair, synthesize, or overwrite scientific records to satisfy this list.
