# Receptor Matching Rules

## 1. Purpose

This document defines the planned matching and grouping rules for comparing structurally similar molecule pairs against public M2OR olfactory-receptor response annotations. The objective is to test whether descriptor-space differences are consistent with receptor-related molecular proxies in public data, not to provide receptor-resolved binding evidence.

## 2. Source database

The receptor annotation analysis uses public M2OR response annotations. Records are interpreted only within the scope of the curated public annotation source and any release-specific filtering documented with the processed tables.

## 3. Molecule matching

- Molecules are matched to M2OR annotations primarily by CAS number.
- When multiple candidate links exist, structure-aware curation is used to avoid identifier mismatch.
- Records with unresolved identity conflicts are excluded from the matched analysis set.

## 4. Duplicate record harmonization

- Duplicate annotations for the same molecule-receptor combination are harmonized according to the curated response annotation available in the public source.
- Response status should be summarized as responsive or non-responsive when such a binary summary is available from the curated record.
- Ambiguous or insufficiently curated duplicates are excluded rather than forced into a definitive category.

## 5. Jointly tested receptor filtering

- For each molecule pair, retain only receptors for which both molecules have curated experimental annotations.
- At least two jointly tested receptors are required for the box-plot analysis and related grouped distance comparisons.
- Pairs with fewer jointly tested receptors may be retained in raw tables for transparency but are excluded from the grouped comparison analysis.

## 6. Pair-level response comparison

- For each retained receptor within a pair, compare the summarized response status of molecule A and molecule B.
- A differential response is counted when the two molecules have different summarized statuses for the same jointly tested receptor.
- The differential-response proportion is calculated as the number of differential responses divided by the number of jointly tested receptors.

## 7. Concordant/divergent classification

- Concordant: no differential responses across the matched receptor subset.
- Divergent: differential-response proportion greater than or equal to 0.5.
- Intermediate: pairs between these two conditions; these are excluded from group comparison.

## 8. Distance comparison

- Euclidean distances are computed separately in standardized QM descriptor space and standardized RDKit descriptor space.
- Standardization parameters must be estimated from the relevant training reference set only.
- Group comparisons are performed only after the jointly tested receptor and concordance/divergence filters described above.

## 9. Interpretation limits

This analysis is an orthogonal consistency analysis. It is designed to assess whether descriptor-space variation is consistent with public receptor-response annotations and receptor-related molecular proxies. It is not receptor-resolved binding evidence, does not prove mechanism, does not reveal receptor coding, and does not directly model receptor combinations.
## Current Local Processing Rules

- Current local input file: `data_pairs_60_pair_OR_comparison.csv`.
- Column names are normalized to snake_case in processed outputs.
- Pair IDs are standardized to `PAIR_0001` style while preserving the original `pair_id_std` in `pair_id_source`.
- CAS matching defaults to normalized CAS columns and falls back to raw CAS when needed.
- Tested status is defined by non-missing responsive annotation after cautious normalization of 1/0/True/False/responsive/nonresponsive strings.
- Duplicate pair-by-species-by-OR records are ranked by jointly tested status, evidence priority, and DOI presence; unresolved conflicts are logged and excluded.
- Concordant pairs have zero differential responses across jointly tested receptors; divergent pairs have differential-response proportion greater than or equal to 0.5; indeterminate pairs are excluded from grouped comparison; fewer than two jointly tested receptors are marked insufficient.
- Human-only outputs are used as the default manuscript-consistent view when Homo sapiens records exist; all-species companion outputs are also retained.
- These public receptor-response annotations are used as an orthogonal consistency analysis and do not establish receptor-resolved binding mechanisms.
