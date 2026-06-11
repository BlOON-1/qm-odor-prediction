# Fig. 2 audit report

## Overall status

- Overall conclusion: PASS with documentation caveats.
- Fig. 2A and Fig. 2B are supported by 24 readable benchmark summary JSON files, and the extracted benchmark AUC values exactly match the raw JSON means.
- Fig. 2C is supported for labelwise mean AUC and delta-AUC calculations across all 24 benchmark combinations, but CNN and MLP JSON files do not expose per-label variance or per-label fold breakdown.
- Fig. 2D is supported by the ten SHAP importance CSV files in `runs_detail`, which were aggregated into `results/shap/shap_feature_ranking.csv`.
- Fig. 2E is supported from the random-forest structure-disjoint summaries in `fig2E`, using the reported mean fields as canonical source values because the repeated `per_fold` entries are placeholders rather than independent folds.

## Release decision

- Promoted into repository: benchmark summary JSON archive, macro-AUC tables, micro-AUC tables, labelwise AUC table, delta-AUC table, SHAP ranking table, structure-disjoint metrics table, and `data/source_data/fig2_source_data.xlsx`.
- Remaining caveats are documentation-level only: CNN/MLP AUC std values were recomputed from benchmark `per_fold` values, and CNN/MLP per-label AUC sources provide means only.
