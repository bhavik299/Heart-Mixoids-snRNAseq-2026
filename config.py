"""
config.py
=========
Centralized configuration, directory paths, publication aesthetics,
and curated gene modules for the Heart Mixoid snRNA-seq pipeline.
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# Base directory defaults to the repository root or an environment override
BASE_DIR = os.environ.get("MIXOID_DATA_DIR", os.path.abspath("./data"))
QC_DIR = os.path.join(BASE_DIR, "QC")
PAIRWISE_DIR = os.path.join(BASE_DIR, "Pairwise_subset")
PUBLIC_DIR = os.path.join(BASE_DIR, "PUBLIC_DATA")
FIG_DIR = os.path.join(BASE_DIR, "FIGURES")
TABLE_DIR = os.path.join(BASE_DIR, "TABLES")

for d in [QC_DIR, PAIRWISE_DIR, FIG_DIR, TABLE_DIR]:
    os.makedirs(d, exist_ok=True)

# Global QC & Clustering Constants
MIN_GENES = 200
MIN_CELLS = 3
MAX_MT_PCT = 5.0       
MAX_GENES_UPPER = 10000
CELL_PRESS_DOUBLE_MM = 174

def mm_to_inches(mm):
    return mm / 25.4

def apply_publication_aesthetics():
    sns.set_style("white")
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "Arial",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.labelweight": "bold",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })
    sc.set_figure_params(dpi=300, facecolor="white", frameon=False)

GREY_REDS_CMAP = LinearSegmentedColormap.from_list(
    "GreyReds", ["#E0E0E0", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"]
)

SAMPLE_PALETTE = {
    "CM14t": "#404040",
    "hSCA06 Normal": "#1F77B4",
    "hSCA06 Hypo": "#D62728",
    "hPSCA02": "#9467BD",
}

ANNOTATION_COLORS = {
    "CM1": "#B2182B", "CM2": "#F46D43", "CM3": "#FFC107",
    "CM4": "#C51B7D", "CM5": "#8E0152",
    "N1": "#08306b", "N2": "#2171b5", "Neuron": "#4292c6", "N3/CM3?": "#9ecae1",
    "Glial/Epicardial": "#33a02c", "Stromal": "#8c564b", "Stromal-Fib": "#e377c2",
}

SUBCLUSTER_PALETTE = [
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#BEBADA",
    "#80B1D3", "#B3DE69", "#FCCDE5", "#FDB462", "#FFFFB3"
]

PAIRWISE_COMPARISONS = {
    "1_hSCA_vs_CM14t": ["hSCA06 Normal", "CM14t"],
    "2_hPSCA_vs_CM14t": ["hPSCA02", "CM14t"],
    "3_Hypoxia_vs_hSCA": ["hSCA06 Normal", "hSCA06 Hypo"],
}

PHYSIOLOGY_MODULES = {
    "Contractile": ["MYH6", "MYH7", "MYL7", "MYL2", "ACTC1"],
    "Structure": ["PDLIM3", "GJA1", "DES"],
    "Immat. Glycolysis": ["LDHA", "ALDOA", "ENO1", "PGK1", "HK1", "HIF1A", "SLC2A1", "PFKP", "PKM"],
    "Immat. OXPHOS": ["COX8A"],
    "Immat. Ion Channels": ["CACNA1H"],
    "Global Maturation": ["ESRRA", "ESRRB", "ESRRG"],
    "Mitochondrial": ["OPA1", "PPARGC1A", "PPARGC1B", "PPARA", "NRF1", "NFE2L2", "MFN1", "MFN2", "TFAM"],
    "Mature Glycolysis": ["SLC2A4", "HK2"],
    "OXPHOS / FAO": ["COX8B", "FABP3", "CPT1B", "ACADVL", "ACADS", "ACADL", "ACADM", "CPT2", "SLC25A20", "HADHA", "HADHB"],
    "Ion Channels": ["SCN5A", "CACNA1C", "CACNB2", "CACNA2D1", "KCNQ1"],
    "Ca2+ (Systole)": ["RYR2"],
    "Ca2+ (Diastole)": ["ATP2A2", "PLN"],
    "Ca2+ (Structure)": ["CAV3", "JPH2", "NEXN", "BIN1", "TRDN"],
}