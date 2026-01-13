# Acupuncture-linked miRNA–mRNA axes map migraine biology and highlight molecular heterogeneity

This repository contains the complete analytical pipeline for the study: **"Acupuncture-linked miRNA–mRNA axes map migraine biology and highlight molecular heterogeneity"**. The workflow integrates literature mining, multi-omics transcriptomics, machine learning, and genetic triangulation to investigate the molecular mechanisms of acupuncture in migraine.

## 1. Overview

The pipeline is organized into five functional modules:
1.  **Network Pharmacology**: Identifying shared targets between acupuncture-responsive mediators and disease-associated genes.
2.  **Transcriptomics**: Differential expression analysis of intervention (GSE198274) and disease-background (E-MTAB-13397, PRJEB40032) cohorts.
3.  **Cross-Omics Integration**: Mapping responsive miRNA–mRNA axes and pathway enrichment.
4.  **Machine Learning**: Subject-aware feature selection and external validation of molecular signatures.
5.  **Causal Inference & Subtyping**: Mendelian randomization and consensus clustering.

## 2. Environment Setup

### Python Environment
```bash
# Create a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### R Environment
This project uses `renv` to manage R dependencies.
```R
# Open R and run:
install.packages("renv")
renv::restore()
```

## 3. Data Acquisition

We provide scripts to automate the download of publicly available datasets.

- **E-MTAB-13397**: Run `python scripts/02_transcriptomics/04_human_emtab13397_acquisition.py` to fetch metadata and count matrices.
- **PRJEB40032**: Run `python scripts/02_transcriptomics/08_prjeb40032_prepare_download_list.py` to generate the ENA download manifest.
- **GSE198274**: Metadata is processed directly via the `05_gse198274_mirna_response.py` script using GEO accessions.

## 4. Analytical Workflow

Execute the scripts in the following order to reproduce the results:

### Phase 1: Target Identification
```bash
python scripts/01_data_mining/01_literature_mining.py
python scripts/01_data_mining/02_target_intersection.py
```

### Phase 2: Transcriptomic Analysis
```R
Rscript scripts/02_transcriptomics/03_meta_dea_by_batch.R
python scripts/02_transcriptomics/05_gse198274_mirna_response.py
python scripts/02_transcriptomics/06_cross_omics_integration.py
```

### Phase 3: Biomarker Discovery
```bash
python scripts/03_ml_biomarkers/02_human_subject_aware_biomarkers.py
python scripts/03_ml_biomarkers/03_external_validation_signature_score.py
```

### Phase 4: Mendelian Randomization & Subtyping
```R
Rscript scripts/04_mr_causality/01_two_sample_mr_opengwas.R
Rscript scripts/05_subtyping/01_consensus_clustering_prjeb40032.R
```

## 5. Directory Structure

- `scripts/`: Source code for all analyses.
- `results/tables/`: Processed statistical tables and performance metrics.
- `figures/final_figures/`: High-resolution figures as presented in the manuscript.
- `docs/`: Supplementary documentation, including image generation prompts for Figure 1.

## 6. License
This project is licensed under the MIT License.
