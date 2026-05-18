# SASA Project 2

## Overview

This repository implements the structural preprocessing pipeline for protein-protein interface analysis in Project 2. The current scope covers two tasks:

- solvent-accessible surface area (`SASA`) computation at atom, residue, and chain levels
- interface residue labeling based on residue-level `ΔSASA`

On top of the core algorithms, the repository also includes a curated dataset of `100` protein complexes, batch-generated interface labels, threshold statistics, and a residue-level training table for downstream machine learning models.

## Objectives

The project uses structural exposure changes to identify interface residues in protein complexes.

The workflow is:

```text
PDB complex
-> SASA calculation
-> apo/holo residue exposure comparison
-> ΔSASA computation
-> binary interface labeling
-> dataset aggregation
-> model-ready residue table
```

In this setting:

- `SASA` measures how much of an atom or residue is accessible to solvent
- `ΔSASA = SASA_apo - SASA_holo`
- residues with sufficiently large `ΔSASA` are treated as interface residues

## Method Summary

### 1. SASA Calculation

The `SASA` module is based on the `Shrake-Rupley` rolling-sphere approximation. It:

- parses atomic coordinates and residue metadata from `PDB`
- loads sphere surface sampling points from `Dot.txt`
- evaluates exposed surface points for each atom
- aggregates atomic surface areas to residue and chain levels

### 2. Interface Label Generation

For a target chain in a protein complex:

- the target chain alone is used to estimate `apo` residue `SASA`
- the target chain together with partner chains is used to estimate `holo` residue `SASA`
- `ΔSASA` is computed residue by residue
- binary labels are generated under multiple thresholds

Thresholds currently used in the repository:

- `0.5`
- `1.0`
- `2.0`
- `5.0`

## Repository Layout

```text
README.md
requirements.txt
data/
  raw/
    examples/                  # example PDB files, dot file, figure assets
    pdb_complexes/             # curated protein complex dataset
  processed/
    examples/                  # example outputs
    interface_labels_per_complex/
    threshold_stats_per_complex/
    complex_manifest.csv
    interface_labels_all.csv
    ml_residue_dataset.csv
    threshold_statistics_by_complex.csv
    threshold_statistics_overall.csv
docs/
  README.md
  module_1_sasa.md
  module_2_delta_sasa.md
  pipeline_overview.md
  project_spec.md
scripts/
  run_sasa_example.py
  run_delta_sasa_example.py
  run_collect_dataset.py
  run_batch_labeling.py
  run_prepare_mlp_dataset.py
src/
  sasa_project/
```

## Core Modules

- [sasa.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/sasa.py:1)
  Implements PDB parsing, dot loading, atom-level SASA computation, and residue/chain aggregation.

- [delta_sasa_label.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/delta_sasa_label.py:1)
  Generates residue-level `ΔSASA` labels for a single complex.

- [batch_generate_interface_labels.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/batch_generate_interface_labels.py:1)
  Applies the labeling pipeline to the full complex dataset and produces aggregate statistics.

- [prepare_mlp_dataset.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/prepare_mlp_dataset.py:1)
  Converts aggregated labeling results into a residue-level table for downstream classification models.

## Data Assets

The repository currently includes:

- `100` curated protein complex structures in [data/raw/pdb_complexes](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/raw/pdb_complexes)
- a complex manifest in [complex_manifest.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/complex_manifest.csv:1)
- residue-level aggregate labels in [interface_labels_all.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_all.csv:1)
- overall threshold statistics in [threshold_statistics_overall.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/threshold_statistics_overall.csv:1)
- a model-ready training table in [ml_residue_dataset.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/ml_residue_dataset.csv:1)

## Usage

Run commands from the repository root with `PYTHONPATH=src`.

### Run the SASA example

```bash
PYTHONPATH=src python3 scripts/run_sasa_example.py
```

### Run the single-complex ΔSASA example

```bash
PYTHONPATH=src python3 scripts/run_delta_sasa_example.py --target-chain C --partner-chains D
```

### Regenerate batch interface labels

```bash
PYTHONPATH=src python3 scripts/run_batch_labeling.py
```

### Regenerate the model-ready residue table

```bash
PYTHONPATH=src python3 scripts/run_prepare_mlp_dataset.py --default-threshold 2.0
```

## Outputs

Key processed outputs are:

- [data/processed/examples](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/examples)
  Example outputs for single-structure SASA and single-complex `ΔSASA`

- [data/processed/interface_labels_per_complex](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_per_complex)
  Per-complex residue labeling results

- [data/processed/threshold_stats_per_complex](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/threshold_stats_per_complex)
  Per-complex threshold statistics

- [interface_labels_all.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_all.csv:1)
  Combined residue-level labels across the dataset

- [ml_residue_dataset.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/ml_residue_dataset.csv:1)
  Residue-level dataset for downstream classification, with the default label set to `ΔSASA > 2.0`

## Documentation

- [docs/README.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/README.md:1)
  Documentation index

- [docs/module_1_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_1_sasa.md:1)
  Notes for the SASA module

- [docs/module_2_delta_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_2_delta_sasa.md:1)
  Notes for the `ΔSASA` labeling module

- [docs/pipeline_overview.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/pipeline_overview.md:1)
  High-level project explanation

## Current Status

The repository currently provides a complete preprocessing pipeline for:

- protein complex collection
- residue-level `apo/holo SASA` computation
- `ΔSASA`-based interface labeling
- threshold-based statistics
- downstream training table preparation

The next natural extension is to attach sequence embeddings, such as `ESM-2`, and train residue-level classifiers on top of the generated structural labels.
