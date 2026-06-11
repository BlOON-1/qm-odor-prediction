# Receptor Annotation Preparation Report

## 1. Input file inspected
- `data_pairs_60_pair_OR_comparison.csv`

## 2. Raw row count
- `3357`

## 3. Pair count
- `60`

## 4. Unique CAS count
- `52`

## 5. Unique OR_gene count
- `849`

## 6. Species distribution
- `Homo sapiens`: `3352`
- `missing`: `5`

## 7. Evidence priority distribution
- `Low`: `2855`
- `Medium`: `432`
- `High`: `56`
- `missing`: `14`

## 8. Jointly tested receptor count distribution by pair
- `0` jointly tested receptors: `28` pairs
- `1` jointly tested receptors: `6` pairs
- `2` jointly tested receptors: `6` pairs
- `3` jointly tested receptors: `1` pairs
- `4` jointly tested receptors: `2` pairs
- `6` jointly tested receptors: `1` pairs
- `7` jointly tested receptors: `1` pairs
- `8` jointly tested receptors: `7` pairs
- `9` jointly tested receptors: `2` pairs
- `10` jointly tested receptors: `4` pairs
- `11` jointly tested receptors: `1` pairs
- `371` jointly tested receptors: `1` pairs

## 9. Number of concordant/divergent/indeterminate/excluded pairs
- concordant: `7`
- divergent: `5`
- indeterminate: `14`
- excluded: `48`

## 10. Duplicate/conflict handling
- conflicting duplicate record groups: `5`
- Evidence priority ranking used: High > Medium > Low > missing.
- Jointly tested records were preferred over untested records.
- DOI-bearing records were preferred over records without DOI when other ranks tied.

## 11. Files generated
- `data/receptor_annotations/m2or_pair_or_raw.csv`
- `data/receptor_annotations/m2or_matched_molecules.csv`
- `data/receptor_annotations/m2or_matched_molecules_all_species.csv`
- `data/receptor_annotations/receptor_response_long.csv`
- `data/receptor_annotations/receptor_response_long_all_species.csv`
- `data/receptor_annotations/receptor_response_vectors.csv`
- `data/receptor_annotations/receptor_response_vectors_all_species.csv`
- `data/receptor_annotations/receptor_pairwise_concordance.csv`
- `data/receptor_annotations/receptor_pairwise_concordance_all_species.csv`
- `data/receptor_annotations/receptor_concordant_divergent_pairs.csv`
- `data/receptor_annotations/receptor_concordant_divergent_pairs_all_species.csv`
- `data/receptor_annotations/receptor_matching_exclusion_log.csv`
- `results/receptor_consistency/receptor_consistency_statistics.csv`

## 12. Files not generated and why
- Descriptor-distance comparison files are handled separately and may be blocked if descriptor matrices are unavailable.

## 13. Interpretation boundary
These public receptor-response annotations are used as an orthogonal consistency analysis and do not establish receptor-resolved binding mechanisms.
