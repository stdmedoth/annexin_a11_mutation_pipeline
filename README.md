# 🧬 Mutation Conformation Pipeline

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/stdmedoth/annexin_a11_mutation_pipeline/blob/main/mutation_pipeline_colab.ipynb)

A computational pipeline for analyzing the conformational effects of protein mutations using BioEmu. This tool generates conformational ensembles for wild-type and mutant proteins, then performs comprehensive structural analysis to understand how mutations affect protein dynamics.

## 🎯 Features

- **Conformational Sampling**: Generate protein conformational ensembles using Microsoft's BioEmu
- **Structural Analysis**: 
  - RMSD/RMSF analysis for flexibility assessment
  - Radius of gyration for compactness evaluation
  - PCA for dimensionality reduction and conformational space visualization
  - Free energy landscape computation
  - Contact map analysis
- **Google Drive Persistence**: Save and load results across Colab sessions
- **Interactive Visualization**: Publication-ready plots for all analyses
- **Binder Design Recommendations**: Automated insights for therapeutic targeting

## 🚀 Quick Start

### Option 1: Google Colab (Recommended)

Click the badge above or go to [Google Colab](https://colab.research.google.com/), then:
1. File → Upload notebook
2. Upload `mutation_pipeline_colab.ipynb`
3. Run All cells

### Option 2: Local Installation

```bash
# Clone the repository
git clone https://github.com/stdmedoth/annexin_a11_mutation_pipeline.git
cd annexin_a11_mutation_pipeline

# Install dependencies
pip install biopython mdtraj prody matplotlib seaborn plotly ipywidgets scipy

# Optional: Install BioEmu for real conformational sampling
pip install git+https://github.com/microsoft/bioemu.git

# Open the notebook
jupyter notebook mutation_pipeline_colab.ipynb
```

## 📋 Workflow

```
1. Input Sequence     →  Enter wild-type protein sequence
2. Define Mutation    →  Specify mutation (e.g., G175E)
3. Generate Conformations  →  BioEmu sampling (or simulated demo data)
4. Structural Analysis     →  RMSD, RMSF, Rg, PCA, Contact Maps
5. Results & Export        →  Figures, summary JSON, downloadable archive
```

## 📊 Analysis Outputs

| Analysis | Description |
|----------|-------------|
| **RMSF** | Per-residue flexibility comparison between WT and mutant |
| **ΔRMSF** | Change in flexibility upon mutation |
| **Radius of Gyration** | Overall protein compactness |
| **PCA** | Conformational space visualization |
| **Free Energy Landscape** | 2D energy surface from PC projections |
| **Contact Maps** | Residue-residue contact probability differences |

## 📁 Output Structure

```
MutationPipeline/
└── ProteinName/
    └── MutationID/
        ├── wt/
        │   ├── *_sequence.fasta
        │   ├── *_input.json
        │   └── *_metadata.json
        ├── mutant/
        │   └── (same structure)
        ├── trajectory_data.npz
        ├── analysis_summary.json
        ├── rmsf_analysis.png
        ├── rg_analysis.png
        ├── pca_analysis.png
        ├── free_energy_landscape.png
        └── contact_maps.png
```

## 🔧 Configuration

### BioEmu Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `num_conformations` | 100 | Number of conformations to generate |
| `device` | cuda | Computing device (cuda/cpu) |
| `seed` | 42 | Random seed for reproducibility |
| `temperature` | 1.0 | Sampling temperature |

## 📖 Example: Annexin A11 G175E

The default example analyzes the G175E mutation in Annexin A11:

```python
protein_name = "AnnexinA11"
mutation_input = "G175E"  # Glycine → Glutamate at position 175
```

## 🔬 Requirements

- Python 3.8+
- numpy
- matplotlib
- seaborn
- scipy
- biopython
- ipywidgets
- BioEmu (optional, for real conformational sampling)

## 📝 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{mutation_pipeline,
  title = {Mutation Conformation Pipeline},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/stdmedoth/annexin_a11_mutation_pipeline}
}
```

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

**Note**: If BioEmu is not installed, the pipeline uses simulated conformational data for demonstration purposes. For production use, install BioEmu following the instructions at [microsoft/bioemu](https://github.com/microsoft/bioemu).
