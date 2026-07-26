# Strain-level Metagenomic Analysis of Butyrate-Producing Burst after FMT in HIV Patients

This repository contains the complete computational pipeline, custom scripts, and source data for the manuscript: 
**"Fecal microbiota transplantation increases intestinal butyrate levels in HIV-infected patients through strain-level functional escalation."** (Replace with your final title if different).

## 📁 Repository Structure

The repository is organized into three main modules corresponding to the Methods section of the paper:

### 1. [Server_Side_Preprocessing_Pipeline](./Server_Side_Preprocessing_Pipeline)
Contains Shell scripts (01-06) for the initial processing of metagenomic data on a Linux server:
- `01_QC_processing.sh`: Raw data quality control using `fastp`.
- `02_assembly_and_binning.sh`: Metagenomic co-assembly and binning using the `MetaWRAP` pipeline.
- `03_taxonomy_and_phylogeny.sh`: Taxonomic identification via `GTDB-Tk` and tree construction using `FastTree`.
- `04_functional_annotation.sh`: Gene prediction and metabolic annotation using `Prokka` and `eggNOG-mapper`.
- `05_genomic_load_calculation.sh`: Calculation of genomic load for butyrate-producing enzymes.
- `06_iRep_calculation.sh`: Calculation of in situ replication rates for core strains.

### 2. [Figure1_iTOL_Source_Files](./Figure1_iTOL_Source_Files)
Contains the raw data used to generate the circular phylogenetic tree in **Figure 1A** via the iTOL platform:
- `Tree_For_iTOL_FINAL.tree`: Phylogenetic tree in Newick format.
- `itol_phylum_strip.txt`, `itol_butyrate_bars.txt`, `itol_hero_symbols.txt`: Annotation files for iTOL visualization.

### 3. [R_Visualization_and_Statistics](./R_Visualization_and_Statistics)
Contains R scripts for statistical analysis and generating Figures 2 through 7:
- `Fig1B_Functional_Centralization_Donut.R`: Donut chart for genomic potential.
- `Fig2_Temporal_Trajectories.R`: Longitudinal abundance tracking.
- `Fig3A_Longitudinal_Boxplot.R`: Total butyrate potential analysis.
- `Fig3B_Pareto_Analysis.R`: Cumulative contribution calculation.
- `Fig5_Genomic_Heatmap.R`: Functional gene profiling.
- `Fig6_Functional_Escalation.R`: Escalation bar chart with step-wise visualization.
- `Fig7_iRep_and_Correlation.R`: iRep comparison and Nadir CD4 clinical correlation.

## 🛠 Dependencies and Requirements

### Bioinformatics Tools
- fastp (v0.21.0)
- Bowtie2 (v2.4.4)
- MetaWRAP (v1.3.2)
- GTDB-Tk (v2.1.1)
- Prokka (v1.14.6)
- eggNOG-mapper (v2.1.9)
- iRep (v1.17)

### R Environment
- R version 4.2.0 or higher
- Required packages: `ggplot2`, `dplyr`, `tidyr`, `readxl`, `pheatmap`, `stringr`.

## 📊 Data Availability
The raw sequencing data analyzed in this study are available in the NCBI SRA repository (BioProject: PRJNA1076935). The 279 high-quality MAGs and processed data tables are archived in Zenodo (DOI: 10.5281/zenodo.21595818).

## ✉️ Contact
For questions regarding the code or the Functional Capacity Weighted Model (FCWM), please contact: [YanLi/1092728032@qq.com]** or [Fengjun Liu/Lfj1162126.com].
