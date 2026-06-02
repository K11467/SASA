# SASA-ESM ICONIP/LNCS Submission Package

This directory contains an anonymous LNCS-style manuscript draft for the SASA project.

## Contents

- `main.tex`: LNCS manuscript source.
- `references.bib`: BibTeX references.
- `llncs.cls`, `splncs04.bst`: Springer LNCS LaTeX class and bibliography style from the CTAN LNCS package.
- `data/experiment_summary_650m.csv`: checkpoint-backed main-corpus metrics.
- `data/leakage_ablation_summary_650m.csv`: strict baselines and holo-aware leakage diagnostics.
- `data/benchmark_dset186_metrics_650m.csv`: Dset_186-local prediction metrics.
- `data/benchmark_pdbtest315_metrics_650m.csv`: PDBtest_315-local prediction metrics.
- `README.md`, `history.txt`: upstream LNCS package notes from CTAN.

## Build

On Overleaf, upload the whole zip and compile `main.tex` with pdfLaTeX.

On a local TeXLive installation:

```bash
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

The source was validated locally with MiKTeX `pdflatex` and `bibtex`.

## Submission Notes

- The manuscript is anonymized. Replace the author/institute block in `main.tex` for camera-ready or non-anonymous submission.
- Reported results are only the experiments completed in this project state.
- The primary model uses ESM-2 plus geometry with no SASA input. Holo-aware rows are explicitly labeled as leakage diagnostics.
- Cross-chain EGNN is retained as an analysis module, not claimed as a stable improvement.
- PDBtest_315-local preprocessing and one-GPU inference are complete: 314 analyzable complexes and 65,119 residues.
- The 1000+ complex and larger ESM model runs are described as planned scale-up work, not as completed results.
