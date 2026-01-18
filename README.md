# Single Cell RNA-seq Analysis Pipeline

A Nextflow-based pipeline for comprehensive single cell RNA-seq analysis using modern tools including Scanpy, SCVI-tools, and Solo for doublet detection.

## Overview

This pipeline performs the following analysis steps:
1. **Data Loading**: Download and format raw count data
2. **AnnData Creation**: Create Scanpy AnnData objects
3. **Doublet Detection**: Remove doublets using Solo
4. **Quality Control**: Filter cells based on QC metrics
5. **Clustering**: Initial clustering with Scanpy
6. **Integration**: Batch correction and integration with SCVI
7. **Differential Expression**: Find marker genes
8. **Cell Type Annotation**: Automated cell type labeling

## Installation

### Prerequisites
- Nextflow (>=22.10.0)
- Docker or Singularity (optional, for containerized execution)
- Conda/Mamba (recommended for dependency management)

### Setup Environment

#### Option 1: Conda Environment (Recommended)
```bash
conda env create -f environment.yml
conda activate single_cell_analysis
```

#### Option 2: pip Environment
```bash
pip install -r requirements.txt
```

#### Option 3: Docker/Singularity
The pipeline uses the `continuumio/miniconda3:latest` container by default, which includes Python and conda for package management.



### ⚠️ IMPORTANT: Correct Command Syntax
**DO NOT USE**: `-params-file samplesheet.csv` ❌
**USE INSTEAD**: `--samplesheet samplesheet.csv` ✅

The `-params-file` option expects JSON/YAML configuration files, not your data samplesheet.

## Input Format

Create a `samplesheet.csv` file with the following format:

```csv
sample_id,counts
sample1,/path/to/sample1_counts.mtx
sample2,/path/to/sample2_counts.mtx
sample3,https://example.com/sample3_counts.mtx
```

**Required columns:**
- `sample_id`: Unique identifier for each sample
- `counts`: Path or URL to count matrix file (MTX, CSV, or H5)

## Output Structure

```
Results/
├── 01_raw_counts/           # Downloaded count files
├── 02_analysis/             # Main analysis results
│   ├── sample1/
│   │   ├── 01_make_adata/   # AnnData objects
│   │   ├── 02_doublets/     # Doublet detection results
│   │   ├── 03_qc_filter/    # QC filtered data
│   │   ├── 04_scanpy_cluster/  # Initial clustering
│   │   ├── 05_scvi_integrate/  # SCVI integration
│   │   ├── 06_scanpy_de_and_plots/  # Differential expression
│   │   └── 07_celltype_annotate_and_plots/  # Cell type annotation
│   └── sample2/
└── pipeline_summary.html    # Execution summary
```

## Configuration

### QC Parameters
Default QC thresholds can be modified by editing the Python scripts:
- **Minimum genes per cell**: 200
- **Maximum genes per cell**: 6000
- **Maximum mitochondrial %**: 20%
- **Maximum ribosomal %**: 5%

### SCVI Parameters
- **Latent dimensions**: 30 (default)
- **Training epochs**: 400 (default)
- **Clustering resolutions**: 0.5 and 1.0

## Resource Requirements

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| DOWNLOAD_COUNTS | 1 | 1 GB | 30m |
| MAKE_ADATA | 2 | 6 GB | 2h |
| SOLO_DOUBLET | 6 | 12 GB | 8h |
| QC_FILTER | 2 | 6 GB | 2h |
| SCANPY_CLUSTER | 4 | 8 GB | 4h |
| SCVI_INTEGRATE | 6 | 14 GB | 10h |
| SCANPY_DE_AND_PLOTS | 4 | 8 GB | 4h |
| CELLTYPE_ANNOTATE_AND_PLOTS | 2 | 4 GB | 2h |



## Advanced Features

### Multi-sample Integration
The pipeline supports batch correction across multiple samples using SCVI. Ensure each sample has a unique identifier in the samplesheet.

### Custom Cell Type Markers
Edit `data/07_celltype_annotate_and_plots.py` to add custom cell type markers for your tissue/cell type.

### Additional Analysis Modules
The modular design allows adding new analysis steps:
1. Create new Python script in `data/`
2. Create corresponding Nextflow module in `Modules/`
3. Include in `main.nf` workflow

## Dependencies

### Python Packages
- scanpy==1.10.1
- scvi-tools>=0.20.0
- anndata>=0.8.0
- pandas>=1.5.0
- numpy>=1.21.0
- matplotlib>=3.5.0
- seaborn>=0.11.0
- scipy>=1.9.0
- scikit-learn>=1.1.0

### System Requirements
- **OS**: Linux, macOS, Windows (with WSL)
- **Python**: 3.9+
- **Nextflow**: 22.10.0+
- **Java**: 11+ (for Nextflow)

