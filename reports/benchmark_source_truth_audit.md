# Benchmark Source Truth Audit

- Source of truth for released benchmark table: `data/processed/benchmark_molecules_4495.csv`.
- Checked benchmark rows: `4495`.
- Raw checked odor vocabulary size: `113`.
- Released manuscript-aligned odor vocabulary size: `112`.
- Release rule used to align the model vocabulary with archived Fig. 2 summary JSON files: exclude `odorless` from the raw 113-label snapshot.
- Conclusion: the current repository now supports a 4,495-molecule benchmark release and a 112-label model vocabulary release, with the label-space reduction explicitly documented.
