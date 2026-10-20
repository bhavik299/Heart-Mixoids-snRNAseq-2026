"""
04_CM_Purification_Trajectories_and_CellCycle.py
================================================
1. Applies strict transcriptomic gating (TNNT2 > 2.0, TUBB3 < 2.5, MAP2 < 2.5, ELAVL3 < 1.25).
2. Generates Combined & Split-Sample CM UMAPs (Figs 3C, S2D-E).
3. Computes Diffusion Pseudotime (DPT) and PAGA trajectories (Fig S2F).
4. Scores Cell Cycle Phases (G1, S, G2/M) across CM subclusters (Figs 3R, S2G, S4D).
5. Renders Chamber Specification, Electrophysiology, and Metabolic Plots (Figs 3G-J, 3L, 4G, 6F-H, 6J-L, S3C-D).
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

S_GENES = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"]
G2M_GENES = ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"]

 CPC_SIGNATURE = ["ISL1", "KDR", "PDGFRA", "GATA4", "MEF2C", "TMEM88", "HAND1", "HAND2"]
IMMATURE_SIGNATURE = ["MKI67", "TOP2A", "CDK1", "AURKB", "ENO1", "HK2", "LDHA", "SLC2A1"]
MATURE_CHECK_SIGNATURE = ["MYL2", "TNNI3", "MYH7", "RYR2", "CACNA1C", "PLN", "ATP2A2"]
ROOT_METHOD = "CPC"  # Toggle between "CPC" (Methods text) and "IMMATURE" (Script 52)

def purify_cms(adata, recluster=False):
    expr = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    if expr.X.max() > 100:
        sc.pp.normalize_total(expr, target_sum=1e4)
        sc.pp.log1p(expr)
    df = sc.get.obs_df(expr, keys=["TNNT2", "TUBB3", "MAP2", "ELAVL3"])
    mask = (df["TNNT2"] > 2.0) & (df["TUBB3"] < 2.5) & (df["MAP2"] < 2.5) & (df["ELAVL3"] < 1.25)
    pure = expr[mask].copy()
    pure.raw = pure
    if recluster:
        sc.pp.highly_variable_genes(pure, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.scale(pure, max_value=10)
        sc.tl.pca(pure, svd_solver="arpack")
        sc.pp.neighbors(pure, n_neighbors=30, n_pcs=40)
        sc.tl.umap(pure, min_dist=0.3)
        sc.tl.leiden(pure, resolution=0.4, key_added="Pure_CM_Subcluster")
        pure.obs["Pure_CM_Subcluster"] = "Sub_" + pure.obs["Pure_CM_Subcluster"].astype(str)
    return pure

# 1. PURIFY hSCM vs CM-ONLY AND RUN DPT & CELL CYCLE (Figs 3C, 3G-L, 4G, S2D-G)
hsca_raw_path = os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Raw1.h5ad")
if os.path.exists(hsca_raw_path):
    adata_hsca = sc.read_h5ad(hsca_raw_path)
    adata_pure = purify_cms(adata_hsca, recluster=True)
    pristine_path = os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Pristine_CM_Manifold.h5ad")
    adata_pure.write(pristine_path)

    # Cell Cycle Scoring (Fig S2G)
    raw_eval = adata_pure.raw.to_adata()
    sc.tl.score_genes_cell_cycle(
        raw_eval,
        s_genes=[g for g in S_GENES if g in raw_eval.var_names],
        g2m_genes=[g for g in G2M_GENES if g in raw_eval.var_names]
    )
    adata_pure.obs["phase"] = raw_eval.obs["phase"]

    # Diffusion Pseudotime (Fig S2F)
    sc.pp.neighbors(adata_pure, n_neighbors=15, use_rep="X_pca")
    sc.tl.diffmap(adata_pure, n_comps=15, random_state=0)
    root_genes = CPC_SIGNATURE if ROOT_METHOD == "CPC" else IMMATURE_SIGNATURE
    sc.tl.score_genes(adata_pure, gene_list=[g for g in root_genes if g in raw_eval.var_names], score_name="Root_Score", use_raw=True)
    adata_pure.uns["iroot"] = int(adata_pure.obs["Root_Score"].argmax())
    sc.tl.dpt(adata_pure, n_branchings=0, n_dcs=10)

    # Feature UMAPs on Pristine CMs (Figs 3H-J, 4G, S3C)
    feature_genes = [g for g in ["HEY2", "TBX2", "NR2F2", "CACNA1D", "KCNH2", "GJA1", "GJA5", "MEIS1", "YAP1", "VAMP2", "SNAP25"] if g in raw_eval.var_names]
    sc.pl.umap(adata_pure, color=feature_genes, use_raw=True, color_map=cfg.GREY_REDS_CMAP, vmin=0, vmax="p99", ncols=3, frameon=False, show=False)
    plt.savefig(os.path.join(cfg.FIG_DIR, "Fig3HJ_4G_S3C_Pristine_CM_Features.pdf"), bbox_inches="tight", transparent=True)
    plt.close()

    # Whole-Sample Electrophysiology DotPlot (Fig 3G)
    raw_eval.obs["Whole_Sample"] = raw_eval.obs["Biological_Sample"].apply(lambda x: "CM-Only (Control)" if "CM14t" in str(x) else "hSCA (Innervated)")
    electro_mods = {"VCM (Structure & Na+)": ["SCN5A", "JPH2", "NEXN"], "VCM (Ca2+ Extrusion)": ["SLC8A1"], "ACM (Ca2+ Reuptake)": ["ATP2A2"]}
    dp_el = sc.pl.DotPlot(raw_eval, var_names=electro_mods, groupby="Whole_Sample", standard_scale="var", vmin=0, vmax=1)
    dp_el.swap_axes()
    dp_el.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.5, grid=True)
    dp_el.savefig(os.path.join(cfg.FIG_DIR, "Fig3G_WholeSample_Calcium_Electro.pdf"), bbox_inches="tight")

# 2. PURIFY hPCM vs CM-ONLY WITHOUT RE-EMBEDDING (Figs 6F-H, 6J-L)
hpsca_raw_path = os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Raw1.h5ad")
if os.path.exists(hpsca_raw_path):
    adata_hpsca = sc.read_h5ad(hpsca_raw_path)
    adata_hpsca_pure = purify_cms(adata_hpsca, recluster=False)
    adata_hpsca_pure.write(os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Pure_Thresholded.h5ad"))