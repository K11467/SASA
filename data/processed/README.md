# Processed Data

This directory contains lightweight tracked summaries and locally generated large artifacts.

## Git Policy

Tracked files:

- benchmark manifests and Delta-SASA labels;
- aggregate metrics and per-residue predictions;
- `experiment_summary_650m.csv` and `leakage_ablation_summary_650m.csv`;
- `artifact_manifest.csv`;
- small example outputs.

Ignored files:

- ESM-2 embedding CSV files;
- multimodal training CSV files;
- downloaded benchmark PDB structures.

The strict no-SASA ESM-2 650M EGNN checkpoint is the primary inference artifact. The
five earlier holo-aware checkpoints are retained for leakage diagnostics and cutoff
analysis. Historical checkpoints remain ignored.

The ignored files are reproducible but too large for normal Git storage. Their expected
sizes and SHA256 values are recorded in `artifact_manifest.csv`.

## PDBtest_315-local Reproduction

```powershell
$env:PYTHONPATH="src"

python scripts\run_benchmark_eval.py --benchmark pdbtest315

python -m sasa_project.extract_esm_embeddings `
  --manifest data\processed\benchmark_pdbtest315_manifest.csv `
  --model-name facebook/esm2_t33_650M_UR50D `
  --device cuda `
  --output data\processed\benchmark_pdbtest315_esm_embeddings_650m.csv

python -m sasa_project.build_multimodal_dataset `
  --labels data\processed\benchmark_pdbtest315_labels.csv `
  --manifest data\processed\benchmark_pdbtest315_manifest.csv `
  --embeddings data\processed\benchmark_pdbtest315_esm_embeddings_650m.csv `
  --output data\processed\benchmark_pdbtest315_multimodal_650m.csv

python -m sasa_project.train_interface_model `
  --input data\processed\benchmark_pdbtest315_multimodal_650m.csv `
  --model egnn `
  --feature-set esm_struct `
  --device cuda `
  --predict-only `
  --checkpoint-input data\processed\best_egnn_650m_esm_struct_d8.pt `
  --prediction-output data\processed\benchmark_pdbtest315_predictions_egnn_esm_struct_650m_d8.csv `
  --metrics-output data\processed\benchmark_pdbtest315_metrics_egnn_esm_struct_650m.csv
```

`PDBtest_315-local` contains 314 analyzable complexes and 65,119 target-chain residues.
Holo-aware EGNN and cross-chain commands remain available for diagnostic analysis, but
their metrics must not be presented as leakage-controlled prediction results.
