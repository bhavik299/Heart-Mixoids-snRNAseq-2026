# Human cardiac organoids innervated by autonomic neurons using a mixoid approach model myocardial infarction under hypoxia

[![DOI](https://zenodo.org/badge/1426878596.svg)](https://doi.org/10.5281/zenodo.23497090)

This repository contains the complete, reproducible Python bioinformatics pipeline for the single-nucleus RNA sequencing (snRNA-seq) analysis reported in **Zeltner et al. (*Cell Stem Cell*, 2026)**.

## 1. System Requirements & Installation
All analyses were executed in Python (v3.10+) using Scanpy (v1.12). To recreate the environment:

```bash
conda create -n mixoids_snrnaseq python=3.10 -y
conda activate mixoids_snrnaseq
pip install -r requirements.txt

2. Pipeline Execution Order & Manuscript Figure Mapping
Scripts are numbered chronologically following the experimental workflow in the STAR Methods.

Script	STAR Methods Section	Manuscript Figures & Tables Generated
01_Primary_QC_and_Global_Integration.py	Primary processing, OCM demultiplexing, QC & global Harmony integration	QC Reproducibility Plots & Global Master Atlas
02_Pairwise_Reintegration_and_Clustree.py	Pairwise manifold reintegration, Clustree/Silhouette optimization & 10k gene cutoff	Clustree & Silhouette Stability Diagnostics
03_Lineage_Annotation_DGE_and_Pathways.py	Unsupervised DGE (Wilcoxon), Enrichr GO:BP 2023 & Lineage Annotation	Figs 3B, 5K–O, 6E, S2A–B, S4A–C, S4E, S5A–B, Table S1
04_CM_Purification_Trajectories_and_CellCycle.py	Strict transcriptomic CM gating (TNNT2 > 2.0), DPT Pseudotime & Cell Cycle	Figs 3C, 3G–J, 3L, 3R, 4G, 6F–H, 6J–L, S2C–G, S2I, S3A, S3C–D, S4D
05_YAP1_MEIS1_Continuum_and_Signaling.py	YAP1/MEIS1 4-state developmental continuum & ERBB2/4 receptor scoring	Figs 3N–Q, 6N–O, S3F–I, S5C–D
06_Reference_Atlas_Benchmarking.py & 42g_Thresholded_Ref_Heatmaps_Spearman.py	Adult Heart Cell Atlas (Pearson) & De Bono Fetal Heart Atlas (Spearman) mapping	Figs 6I, 6P, S2H, S4F

3. Note on Figure Assembly
The scripts in this repository execute all statistical modeling, gating, and base .pdf figure generation. For the final publication, base vector PDFs were imported into Adobe Illustrator solely for multi-panel layout alignment and typography standardization without altering underlying data coordinates.

4. Data & Code Availability
Raw & Processed Sequencing Data: Deposited at NCBI GEO under accession number GSE341502.

Lead Contact: Further information and requests should be directed to and will be fulfilled by the Lead Contact, Dr. Nadja Zeltner (nadja.zeltner@uga.edu).