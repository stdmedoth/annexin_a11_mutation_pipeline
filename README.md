# Mutation Conformation Pipeline

A pipeline for generating and analyzing protein conformations for mutations using BioEmu.

## Overview

This pipeline allows you to:
1. Input a wild-type (WT) protein sequence
2. Specify mutations to study
3. Generate conformational ensembles using BioEmu
4. Align and analyze the conformations

## Project Structure

```
annexin_a11_mutation_pipeline/
├── config/
│   └── settings.yaml          # Configuration settings
├── src/
│   ├── __init__.py
│   ├── sequence_handler.py    # Sequence manipulation utilities
│   ├── mutation_handler.py    # Mutation parsing and application
│   ├── bioemu_runner.py       # BioEmu integration
│   ├── alignment.py           # Structure alignment utilities
│   └── analysis.py            # Conformational analysis tools
├── data/
│   ├── sequences/             # Input WT sequences
│   ├── mutations/             # Mutation definitions
│   └── outputs/               # Generated conformations
├── scripts/
│   └── run_pipeline.py        # Main pipeline script
├── requirements.txt
└── README.md
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python scripts/run_pipeline.py
```

## Analysis Methods

After generating conformations, you can analyze:
- RMSD/RMSF profiles
- Radius of gyration
- Secondary structure content
- Contact maps
- Free energy landscapes
- Principal Component Analysis (PCA)
