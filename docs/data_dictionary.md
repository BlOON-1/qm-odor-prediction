# Data Dictionary

This document summarizes the manuscript-aligned release tables together with the currently checked raw local snapshot files.

## Manuscript-aligned release files

| Target file | Purpose | Key columns or structure | Notes |
| --- | --- | --- | --- |
| `benchmark_molecules_4495.csv` | Processed benchmark molecule table for 4,495 curated odorants. | `molecule_id`, `cas_number`, `canonical_smiles`, `source_record_group`, `structure_verification_status`, `included_in_benchmark` | Molecules are aligned by identifier and verified against structure where possible; unresolved conflicts are excluded. |
| `odor_label_vocabulary_112.csv` | Standardized benchmark odor vocabulary used for modeling. | `label_id`, `standard_label`, `normalized_label`, `label_group`, `label_definition_note` | This release vocabulary excludes the `odorless` label from the raw 113-label snapshot so that the public table matches the model-training vocabulary used in Fig. 2 benchmark summaries. |
| `odor_label_matrix_4495x112.csv` | Multi-hot odor label matrix for the benchmark set. | `molecule_id` plus 112 binary label columns | Label columns correspond exactly to the released 112-label vocabulary table. |
| `label_harmonization_table.csv` | Mapping from source labels to standardized labels. | `source_name`, `raw_label`, `normalized_label`, `standard_label`, `mapping_rule`, `notes` | Documents vocabulary normalization and harmonization decisions. |
| `rdkit_descriptors.csv` | RDKit-computed 2D descriptor matrix. | `molecule_id`, descriptor columns | Continuous features are standardized within downstream training workflows, not in the raw release table. |
| `ecfp4_fingerprints.npz` | Binary or sparse ECFP4 fingerprint storage. | Array-like fingerprint matrix plus `molecule_id` index mapping | Stored as a compact matrix rather than CSV because of dimensionality. |
| `iterative_stratified_splits/*.csv` | Train-test indices for repeated iterative multi-label stratification. | `split_name`, `seed`, `molecule_id`, `partition` | Planned seeds are 42 through 46 with a 7:3 train-test split. |
| `structure_disjoint_splits/*.csv` | Structure-disjoint split assignments. | `split_name`, `seed`, `molecule_id`, `partition`, `max_train_similarity` | High-similarity structural neighbors are defined by ECFP4 Tanimoto similarity greater than 0.8. |
| `external_molecules_87.csv` | External fragrance molecule metadata for sensory validation. | `external_id`, `molecule_id_or_alias`, `canonical_smiles`, `pair_group`, `availability_notes` | Covers 87 external molecules used for the manuscript-aligned external comparison set. |
| `structural_neighbor_pairs_60.csv` | Structural-neighbor pair list used for external comparisons. | `pair_id`, `molecule_a`, `molecule_b`, `pair_selection_rule`, `similarity_metric`, `similarity_value` | Release file for 60 near-neighbor pairs. |
| `panel_label_vocabulary_43.csv` | Sensory-panel vocabulary derived from the benchmark vocabulary. | `panel_label_id`, `panel_label`, `benchmark_label_link`, `notes` | Contains the 43-label panel vocabulary. |
| `aggregated_panel_ratings.csv` | Aggregated sensory-panel response summaries. | `molecule_id`, `panel_label`, `mean_intensity`, `selection_frequency`, `n_assessors` | Only anonymized and aggregated summaries are released. |
| `panel_derived_topk_labels.csv` | Dominant label sets derived from ordered mean-intensity profiles. | `molecule_id`, `rank_order`, `panel_label`, `mean_intensity`, `knee_point_included` | Top-K labels are selected until a knee point in the mean-intensity curve. |
| `pairwise_panel_jaccard_similarity.csv` | Pairwise perceptual similarity between external molecule pairs. | `pair_id`, `molecule_a`, `molecule_b`, `jaccard_similarity`, `label_set_size_a`, `label_set_size_b` | Similarity is defined on panel-derived label sets. |
| `receptor_response_vectors.csv` | Matched public receptor-response annotation vectors. | `pair_id`, `molecule_id`, `receptor_id`, `response_status`, `curation_note` | Includes only matched and curated public annotations used in the orthogonal consistency analysis. |
| `receptor_pairwise_concordance.csv` | Pair-level receptor comparison summaries. | `pair_id`, `n_jointly_tested_receptors`, `n_differential_responses`, `differential_response_proportion`, `classification` | Intermediate pairs may be retained here even if excluded from grouped comparisons. |
| `receptor_concordant_divergent_pairs.csv` | Final receptor concordance and divergence grouping table. | `pair_id`, `classification`, `included_in_group_comparison` | Intended for group-level comparisons after filtering. |

## General notes

- Stable molecule identifiers should be consistent across all released tables.
- The raw checked vocabulary snapshot presently contains 113 labels, while the manuscript-aligned release vocabulary contains 112 labels because `odorless` is excluded from the model-training target space.
- Exact column order may vary slightly across final files, but field meaning will remain consistent with this dictionary.

## Currently generated local snapshot files

| File | Key fields | Purpose |
| --- | --- | --- |
| `data/processed/benchmark_molecules_observed.csv` | `molecule_id`, `cas`, `canonical_smiles`, `odor_labels_cleaned` | Current checked raw benchmark snapshot with 4,495 rows. |
| `data/processed/benchmark_molecules_4495.csv` | `molecule_id`, `cas`, `canonical_smiles`, `odor_labels_cleaned` | Manuscript-aligned benchmark release table with 4,495 rows. |
| `data/processed/odor_label_vocabulary_observed.csv` | `label_id`, `raw_label`, `normalized_label`, `label_frequency` | Current checked raw odor-vocabulary snapshot with 113 rows. |
| `data/processed/odor_label_vocabulary_112.csv` | `label_id`, `raw_label`, `normalized_label`, `label_frequency` | Manuscript-aligned model vocabulary with `odorless` excluded. |
| `data/processed/odor_label_matrix_observed.csv` | `molecule_id` plus label columns | Current checked raw multi-hot benchmark label matrix with 113 label columns. |
| `data/processed/odor_label_matrix_4495x112.csv` | `molecule_id` plus 112 binary label columns | Manuscript-aligned benchmark label matrix matching the released 112-label model vocabulary. |
| `data/external_panel/external_molecules_87.csv` | `external_molecule_id`, `cas`, `molecule_id`, `canonical_smiles` | External molecule metadata table matching the manuscript-aligned 87-molecule set. |
| `data/external_panel/aggregated_panel_ratings.csv` | `cas`, `label`, summary statistics | Aggregated sensory-panel ratings by molecule and label. |
| `data/external_panel/panel_derived_topk_labels.csv` | `cas`, `top1_label`-`top5_label`, scores | Dominant labels derived from panel mean intensities. |
| `data/external_panel/structural_neighbor_pairs_60.csv` | `pair_id`, `cas1`, `cas2`, `jaccard_similarity`, `tanimoto_similarity` | Structural-neighbor pair table derived from the current local input. |
| `data/source_data/fig1_source_data.xlsx` | multiple sheets | Figure 1 source data workbook built from the current checked benchmark snapshot. |
| `data/source_data/fig3_source_data.xlsx` | multiple sheets | Figure 3 source data workbook built from the current checked panel and pair snapshot. |
| `data/source_data/fig4_source_data.xlsx` | multiple sheets | Figure 4 source data workbook built from the current checked receptor-analysis snapshot. |

## RDKit and ECFP4 representation files

| File | Key fields or structure | Purpose |
| --- | --- | --- |
| `data/processed/rdkit_descriptors_raw.csv` | `molecule_id`, `cas`, `chemical_name`, `canonical_smiles`, RDKit descriptor columns | Raw RDKit 2D continuous descriptor table after removal of descriptor columns that are entirely missing. |
| `data/processed/rdkit_descriptors_unstandardized_qc.csv` | identifier columns plus RDKit descriptors passing global missingness and constant-column QC | QC-only RDKit descriptor table produced before any train/test-aware standardization. |
| `data/processed/ecfp4_fingerprints.npz` | sparse binary matrix | Compressed 1024-bit Morgan fingerprint matrix with radius 2. |
| `data/processed/ecfp4_fingerprint_index.csv` | `molecule_id`, `cas`, `chemical_name`, `canonical_smiles`, `row_index` | Row mapping between molecules and the ECFP4 sparse matrix. |
| `data/processed/ecfp4_fingerprint_metadata.json` | settings and molecule counts | Audit metadata for ECFP4 generation parameters and validity counts. |

## Receptor and structure-disjoint files

| File | Key fields or structure | Purpose |
| --- | --- | --- |
| `data/receptor_annotations/m2or_pair_or_raw.csv` | raw M2OR pairwise comparison rows | Local clean copy of the uploaded receptor-comparison source table. |
| `data/receptor_annotations/m2or_matched_molecules.csv` | `cas`, `molecule_role_source`, count summaries | Molecule-level summary of matched receptor annotation coverage under the default species filter. |
| `data/receptor_annotations/receptor_response_long.csv` | `pair_id`, `or_gene`, tested flags, response flags | Harmonized pair-receptor annotation table used for final comparison logic. |
| `data/receptor_annotations/receptor_response_vectors.csv` | one row per CAS with OR-gene columns | Molecule-by-receptor response matrix for the default species filter. |
| `data/receptor_annotations/receptor_pairwise_concordance.csv` | pair-level jointly tested receptor summary columns | Pairwise receptor concordance table used for orthogonal consistency grouping. |
| `data/receptor_annotations/receptor_concordant_divergent_pairs.csv` | concordant and divergent pairs only | Final grouped pair table retained for group comparison. |
| `data/receptor_annotations/receptor_matching_exclusion_log.csv` | exclusion reason log | Explicit log of rows excluded during duplicate handling, human-only filtering, or missing-data filtering. |
| `data/splits/structure_disjoint_splits/structure_disjoint_train.csv` | split table with `cas` | Canonical structure-disjoint training split organized from uploaded fold files. |
| `data/splits/structure_disjoint_splits/structure_disjoint_test.csv` | split table with `cas` | Canonical structure-disjoint test split organized from uploaded fold files. |
| `data/splits/structure_disjoint_splits/structure_disjoint_split_audit.csv` | audit metrics | Leakage and role-convention audit for the uploaded structure-disjoint split files. |
