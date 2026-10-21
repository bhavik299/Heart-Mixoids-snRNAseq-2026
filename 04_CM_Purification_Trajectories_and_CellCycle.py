"""
04_CM_Purification_Trajectories_and_CellCycle.py
================================================
Complete Cardiomyocyte Purification, Trajectory, Electrophysiology, and RYR2 Suite.

Generates:
1. Strict transcriptomic CM gating (TNNT2 > 2.0, TUBB3 < 2.5, MAP2 < 2.5, ELAVL3 < 1.25).
2. Combined & Background-Shadow Split CM UMAPs (Figs 3C, S2D-E).
3. Diffusion Pseudotime (DPT), PAGA, & Cell Cycle Phase Bars (Figs 3R, S2F, S2G, S4D).
4. Mature Sub_3 vs Sub_2 & Trajectory Electrophysiology DotPlots (Fig 3G).
5. Chamber TFs, Synaptic Features, & Metabolic Heatmaps (Figs 3H-J, 3L, 4G, S3C, S3D).
6. hPCM vs CM-Only Thresholded UMAPs (TNNT2, RYR2, ACM/VCM), RYR2+ Percentages,
   and Targeted Maturation Heatmaps/DotPlots (Figs 6F, 6G, 6H, 6J, 6K, 6L).
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

# Standard Cell Cycle Genes (Tirosh et al., 2016)
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8"
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
]

CPC_SIGNATURE = ["ISL1", "KDR", "PDGFRA", "GATA4", "MEF2C", "TMEM88", "HAND1", "HAND2"]
IMMATURE_SIGNATURE = ["MKI67", "TOP2A", "CDK1", "AURKB", "ENO1", "HK2", "LDHA", "SLC2A1"]
ROOT_METHOD = "CPC"


def purify_cms(adata, recluster=False):
    """Applies strict in silico CM gating (TNNT2 > 2.0, neural markers excluded)."""
    obsm_backup = {k: adata.obsm[k].copy() for k in adata.obsm.keys()}
    expr = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    for k, v in obsm_backup.items():
        expr.obsm[k] = v

    if expr.X.max() > 100:
        sc.pp.normalize_total(expr, target_sum=1e4)
        sc.pp.log1p(expr)

    df = sc.get.obs_df(expr, keys=["TNNT2", "TUBB3", "MAP2", "ELAVL3"])
    mask = (
        (df["TNNT2"] > 2.0) &
        (df["TUBB3"] < 2.5) &
        (df["MAP2"] < 2.5) &
        (df["ELAVL3"] < 1.25)
    )
    pure = expr[mask].copy()
    pure.raw = pure

    if recluster:
        pure.obsm.clear()
        pure.varm.clear()
        sc.pp.highly_variable_genes(pure, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.scale(pure, max_value=10)
        sc.tl.pca(pure, svd_solver="arpack")
        sc.pp.neighbors(pure, n_neighbors=30, n_pcs=40)
        sc.tl.umap(pure, min_dist=0.3)
        sc.tl.leiden(pure, resolution=0.4, key_added="Pure_CM_Subcluster")
        pure.obs["Pure_CM_Subcluster"] = "Sub_" + pure.obs["Pure_CM_Subcluster"].astype(str)
    return pure


# =============================================================================
# PART 1: hSCM vs CM-ONLY PRISTINE MANIFOLD (Figs 3C, 3G-L, 4G, S2D-G, S3C-D)
# =============================================================================
hsca_raw_path = os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Raw1.h5ad")
if os.path.exists(hsca_raw_path):
    adata_hsca = sc.read_h5ad(hsca_raw_path)
    adata_pure = purify_cms(adata_hsca, recluster=True)
    pristine_path = os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Pristine_CM_Manifold.h5ad")
    adata_pure.write(pristine_path)

    raw_eval = adata_pure.raw.to_adata()
    raw_eval.obsm["X_umap"] = adata_pure.obsm["X_umap"].copy()

    # 1A. Combined & Split Background-Shadow UMAPs (Figs 3C, S2D-E)
    umap_coords = adata_pure.obsm["X_umap"]
    unique_samples = sorted(adata_pure.obs["Biological_Sample"].dropna().unique().tolist())
    fig_split, axes_split = plt.subplots(1, len(unique_samples), figsize=(6 * len(unique_samples), 6))
    if len(unique_samples) == 1:
        axes_split = [axes_split]

    for idx, s_name in enumerate(unique_samples):
        ax = axes_split[idx]
        ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c="#E0E0E0", s=1.5, alpha=0.6, edgecolors="none", zorder=0)
        s_mask = (adata_pure.obs["Biological_Sample"] == s_name).values
        sc.pl.umap(adata_pure[s_mask], color="Pure_CM_Subcluster", ax=ax, frameon=False, show=False, title=s_name)
    fig_split.savefig(os.path.join(cfg.FIG_DIR, "Fig3C_S2DE_Pristine_Split_UMAP.pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig_split)

    # 1B. Cell Cycle Scoring & Diffusion Pseudotime (Figs S2F, S2G)
    sc.tl.score_genes_cell_cycle(
        raw_eval,
        s_genes=[g for g in S_GENES if g in raw_eval.var_names],
        g2m_genes=[g for g in G2M_GENES if g in raw_eval.var_names],
    )
    adata_pure.obs["phase"] = raw_eval.obs["phase"]

    sc.pp.neighbors(adata_pure, n_neighbors=15, use_rep="X_pca")
    sc.tl.diffmap(adata_pure, n_comps=15, random_state=0)
    root_genes = CPC_SIGNATURE if ROOT_METHOD == "CPC" else IMMATURE_SIGNATURE
    sc.tl.score_genes(adata_pure, gene_list=[g for g in root_genes if g in raw_eval.var_names], score_name="Root_Score", use_raw=True)
    adata_pure.uns["iroot"] = int(adata_pure.obs["Root_Score"].argmax())
    sc.tl.dpt(adata_pure, n_branchings=0, n_dcs=10)

    # 1C. Feature UMAPs on Pristine CMs (Figs 3H-J, 4G, S3C)
    feature_genes = [
        g for g in ["HEY2", "TBX2", "NR2F2", "CACNA1D", "KCNH2", "GJA1", "GJA5", "MEIS1", "YAP1", "VAMP2", "SNAP25"]
        if g in raw_eval.var_names
    ]
    sc.pl.umap(adata_pure, color=feature_genes, use_raw=True, color_map=cfg.GREY_REDS_CMAP, vmin=0, vmax="p99", ncols=3, frameon=False, show=False)
    plt.savefig(os.path.join(cfg.FIG_DIR, "Fig3HJ_4G_S3C_Pristine_CM_Features.pdf"), bbox_inches="tight", transparent=True)
    plt.close()

    # 1D. Mature Sub_3 vs Sub_2 Electrophysiology & TF DotPlots (Figs 3G, S3D)
    def map_mature_subclusters(row):
        if row["Biological_Sample"] == "CM14t" and row["Pure_CM_Subcluster"] == "Sub_3":
            return "CM14t (mature Sub_3)"
        if row["Biological_Sample"] == "hSCA06 Normal" and row["Pure_CM_Subcluster"] == "Sub_2":
            return "hSCA (mature Sub_2)"
        return "Exclude"

    raw_eval.obs["Mature_Subcluster"] = raw_eval.obs.apply(map_mature_subclusters, axis=1)
    sub_mature = raw_eval[raw_eval.obs["Mature_Subcluster"] != "Exclude"].copy()
    sub_mature.obs["Mature_Subcluster"] = pd.Categorical(
        sub_mature.obs["Mature_Subcluster"],
        categories=["CM14t (mature Sub_3)", "hSCA (mature Sub_2)"],
        ordered=True,
    )

    # Fig 3G: 5-Gene Core Calcium & Electrophysiology Panel (Sub_3 vs Sub_2)
    electro_5gene = {
        "VCM (Structure & Na+)": ["SCN5A", "JPH2", "NEXN"],
        "VCM (Ca2+ Extrusion)": ["SLC8A1"],
        "ACM (Ca2+ Reuptake)": ["ATP2A2"],
    }
    dp_el5 = sc.pl.DotPlot(sub_mature, var_names=electro_5gene, groupby="Mature_Subcluster", standard_scale="var", vmin=0, vmax=1, figsize=(3.2, 4.0))
    dp_el5.swap_axes()
    dp_el5.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.8, grid=True)
    dp_el5.savefig(os.path.join(cfg.FIG_DIR, "Fig3G_Sub3_vs_Sub2_5Gene_Calcium_Electro.pdf"), bbox_inches="tight")

    # Fig 3G (Extended 9-Gene VCM vs ACM Panel on Sub_3 vs Sub_2)
    electro_9gene = {
        "Ventricular Features (VCM)": ["SCN5A", "KCNH2", "JPH2", "NEXN", "SLC8A1", "GJA1"],
        "Atrial Features (ACM)": ["CACNA1D", "ATP2A2", "GJA5"],
    }
    dp_el9 = sc.pl.DotPlot(sub_mature, var_names=electro_9gene, groupby="Mature_Subcluster", standard_scale="var", vmin=0, vmax=1, figsize=(5.5, 2.6))
    dp_el9.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.8, grid=True)
    dp_el9.savefig(os.path.join(cfg.FIG_DIR, "Fig3G_Sub3_vs_Sub2_9Gene_VCM_vs_ACM.pdf"), bbox_inches="tight")

    # Fig S3D: Targeted Atrial vs Ventricular TFs (Sub_3 vs Sub_2)
    tf_mods = {"Atrial TFs": ["HEY1", "BMP2", "TBX5"], "Ventricular TFs": ["IRX4", "HAND1", "IRX3"]}
    dp_tf = sc.pl.DotPlot(sub_mature, var_names=tf_mods, groupby="Mature_Subcluster", standard_scale="var", vmin=0, vmax=1, figsize=(5.0, 2.5))
    dp_tf.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=1.0, grid=True)
    dp_tf.savefig(os.path.join(cfg.FIG_DIR, "FigS3D_Targeted_ACM_vs_VCM_TFs.pdf"), bbox_inches="tight")

    # Fig 3L: Head-to-Head Mature Metabolic Heatmap (Sub_3 vs Sub_2)
    metab_mods = {
        "Global Maturation": ["ESRRA", "ESRRB", "ESRRG"],
        "Mitochondrial": ["OPA1", "PPARGC1A", "PPARGC1B", "PPARA", "NRF1", "NFE2L2", "MFN1", "MFN2", "TFAM"],
        "OXPHOS / FAO": ["COX8B", "FABP3", "CPT1B", "ACADVL", "ACADS", "ACADL", "ACADM", "CPT2", "SLC25A20", "HADHA", "HADHB"],
    }
    valid_metab = {k: [g for g in v if g in sub_mature.var_names] for k, v in metab_mods.items()}
    mp_metab = sc.pl.MatrixPlot(sub_mature, valid_metab, groupby="Mature_Subcluster", standard_scale="var", figsize=(3.5, 6.0))
    mp_metab.swap_axes()
    mp_metab.style(cmap="viridis", edge_color="white", edge_lw=0.5)
    mp_metab.savefig(os.path.join(cfg.FIG_DIR, "Fig3L_Mature_Metabolic_MatrixPlot.pdf"), bbox_inches="tight")


# =============================================================================
# PART 2: hPCM vs CM-ONLY THRESHOLDED CMs & RYR2 SUITE (Figs 6F-H, 6J-L)
# =============================================================================
hpsca_raw_path = os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Raw1.h5ad")
if os.path.exists(hpsca_raw_path):
    adata_hpsca = sc.read_h5ad(hpsca_raw_path)
    adata_hpsca_pure = purify_cms(adata_hpsca, recluster=False)
    thresh_path = os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Pure_Thresholded.h5ad")
    adata_hpsca_pure.write(thresh_path)

    # 2A. Masked Feature UMAPs for TNNT2, RYR2, and ACM/VCM Markers (Figs 6F, 6G, 6K)
    adata_plot_u = adata_hpsca.raw.to_adata() if adata_hpsca.raw is not None else adata_hpsca.copy()
    for k in adata_hpsca.obsm.keys():
        adata_plot_u.obsm[k] = adata_hpsca.obsm[k].copy()
    if "X_umap_15" in adata_plot_u.obsm:
        adata_plot_u.obsm["X_umap"] = adata_plot_u.obsm["X_umap_15"].copy()

    adata_plot_u = adata_plot_u[adata_plot_u.obs["Manuscript_Identity"].isin(["CM1", "CM2", "CM4"])].copy()
    if adata_plot_u.X.max() > 100:
        sc.pp.normalize_total(adata_plot_u, target_sum=1e4)
        sc.pp.log1p(adata_plot_u)

    tnnt2_vec = sc.get.obs_df(adata_plot_u, keys=["TNNT2"])["TNNT2"].values
    valid_cm_mask = tnnt2_vec > 2.0

    fig6_markers = ["TNNT2", "RYR2", "MYL2", "MYL7", "CACNA1C", "ATP2A2", "NR2F2", "CACNA1D", "HEY2", "KCNH2"]
    for gene in fig6_markers:
        if gene not in adata_plot_u.var_names:
            continue
        g_vec = sc.get.obs_df(adata_plot_u, keys=[gene])[gene].values
        adata_plot_u.obs[f"masked_{gene}"] = np.where(valid_cm_mask, g_vec, np.nan)

        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.umap(
            adata_plot_u, color=f"masked_{gene}", color_map=cfg.GREY_REDS_CMAP,
            na_color="#E0E0E0", vmin=0, vmax="p99", sort_order=True,
            frameon=False, title=f"{gene} (TNNT2 > 2)", ax=ax, show=False
        )
        fig.savefig(os.path.join(cfg.FIG_DIR, f"Fig6FGK_Thresholded_UMAP_{gene}.pdf"), dpi=300, bbox_inches="tight", transparent=True)
        plt.close(fig)
        del adata_plot_u.obs[f"masked_{gene}"]

    # 2B. RYR2+ Cardiomyocyte Percentage Calculation & Barplot (Fig 6H)
    adata_cm_thresh = adata_plot_u[valid_cm_mask].copy()
    ryr2_pos = sc.get.obs_df(adata_cm_thresh, keys=["RYR2"])["RYR2"].values > 0
    adata_cm_thresh.obs["RYR2_Positive"] = ryr2_pos
    adata_cm_thresh.obs["Line"] = adata_cm_thresh.obs["Biological_Sample"].astype(str).apply(
        lambda x: "hPCM" if "hPSCA" in x else "CM-only"
    )
    ryr2_pct = adata_cm_thresh.obs.groupby("Line")["RYR2_Positive"].mean() * 100
    ryr2_pct.to_csv(os.path.join(cfg.TABLE_DIR, "Fig6H_RYR2_Positive_Percentages.csv"))

    fig, ax = plt.subplots(figsize=(3.5, 4.5))
    ryr2_pct.reindex(["CM-only", "hPCM"]).plot(kind="bar", color=["#404040", "#9467BD"], edgecolor="black", width=0.6, ax=ax)
    ax.set_ylabel("% of CMs Expressing RYR2", fontweight="bold")
    ax.set_ylim(0, 100)
    sns.despine()
    fig.savefig(os.path.join(cfg.FIG_DIR, "Fig6H_RYR2_Percentage_Barplot.pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig)

    # 2C. Targeted CM Maturation Heatmap & DotPlot (Figs 6J, 6L: CM14t-CM1, CM14t-CM2, hPSCA-CM4)
    adata_cm_thresh.obs["Condition"] = (
        adata_cm_thresh.obs["Biological_Sample"].astype(str).apply(lambda x: "hPSCA" if "hPSCA" in x else "CM14t")
        + " - " + adata_cm_thresh.obs["Manuscript_Identity"].astype(str)
    )
    target_conds = ["CM14t - CM1", "CM14t - CM2", "hPSCA - CM4"]
    sub_6jl = adata_cm_thresh[adata_cm_thresh.obs["Condition"].isin(target_conds)].copy()
    sub_6jl.obs["Condition"] = pd.Categorical(sub_6jl.obs["Condition"], categories=target_conds, ordered=True)

    maturation_dict = {
        "Master Metabolism": ["ESRRA", "ESRRB", "ESRRG", "YAP1", "MEIS1"],
        "Mitochondria": [
            "OPA1", "PPARGC1A", "PPARGC1B", "PPARA", "NRF1", "MFN1", "MFN2",
            "NFE2L2", "TFAM", "CPT1B", "ACADVL", "ACADS", "ACADL",
            "ACADM", "CPT2", "SLC25A20", "HADHA", "HADHB"
        ],
        "Mature Glycolysis": ["SLC2A4", "HK2"],
        "Ion Channels": ["SCN5A", "CACNA1C", "CACNB2", "CACNA2D1", "KCNQ1"],
        "Calcium Systole/Diastole": ["RYR2", "ATP2A2", "PLN"],
        "Calcium Structure": ["CAV3", "JPH2", "NEXN", "BIN1", "TRDN"],
    }
    clean_mat_dict = {k: [g for g in v if g in sub_6jl.var_names] for k, v in maturation_dict.items() if any(g in sub_6jl.var_names for g in v)}
    dyn_h = max(8.0, sum(len(v) for v in clean_mat_dict.values()) * 0.22)

    mp_6l = sc.pl.matrixplot(sub_6jl, var_names=clean_mat_dict, groupby="Condition", standard_scale="var", swap_axes=True, figsize=(4, dyn_h), return_fig=True)
    mp_6l.style(edge_color="black", edge_lw=0.5)
    mp_6l.savefig(os.path.join(cfg.FIG_DIR, "Fig6L_CM_Maturation_Targeted_Heatmap.pdf"), bbox_inches="tight")

    dp_6j = sc.pl.dotplot(sub_6jl, var_names=clean_mat_dict, groupby="Condition", standard_scale="var", swap_axes=True, figsize=(4, dyn_h), return_fig=True)
    dp_6j.style(cmap="viridis", dot_edge_color="black", dot_edge_lw=0.5)
    dp_6j.savefig(os.path.join(cfg.FIG_DIR, "Fig6J_CM_Maturation_Targeted_DotPlot.pdf"), bbox_inches="tight")
