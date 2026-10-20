"""
03_Lineage_Annotation_DGE_and_Pathways.py
=========================================
Generates:
- Annotated Global & Split-Sample UMAPs + Cluster Composition Bars (Figs 3B, 5K, 5L, 6E, S2B, S4B, S5A)
- Master 16-Lineage Annotation DotPlots (Figs S2A, S4A, S5B)
- Balanced Wilcoxon DGE Marker Tables (Table S1)
- Per-Cluster & Butterfly GO Biological Process 2023 Enrichments (Figs 5O, S4C)
- Hypoxia Apoptosis & NRG1/ERBB4 MatrixPlots (Figs 5M, 5N, S4E)
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gseapy as gp
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

PAIRWISE_FILES = {
    "1_hSCA_vs_CM14t": os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Raw1.h5ad"),
    "2_hPSCA_vs_CM14t": os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Raw1.h5ad"),
    "3_Hypoxia_vs_hSCA": os.path.join(cfg.PAIRWISE_DIR, "Hypoxia_vs_hSCA_Raw1.h5ad"),
}

LINEAGE_MARKERS = {
    "Pan-Stromal": ["PDGFRA"],
    "Fibroblast": ["POSTN", "COL1A1"],
    "Epicardial": ["WT1", "TBX18", "ALDH1A2"],
    "Endothelial": ["PECAM1"],
    "SMC": ["ACTA2", "TAGLN", "CNN1"],
    "Neural Progenitor": ["NES", "SOX2", "ASCL1", "GATA2", "GATA3"],
    "Pan-Neural": ["CACNA1B", "SNAP25", "STMN2", "VAMP2", "MAP2", "TUBB3", "MEIS1", "SYT1", "NRXN1", "CHGA"],
    "Glial/Schwann": ["S100B", "MBP", "ERBB3", "PLP1", "FABP7", "ERBB4"],
    "Developing PNS": ["PHOX2A", "PHOX2B", "GFRA3"],
    "Sympathetic": ["TH"],
    "Parasympathetic": ["CHAT", "SLC5A7"],
    "Peripheral Neuron": ["PRPH"],
    "Neural Receptors": ["NTRK1", "NTRK2", "NTRK3", "GFRA1", "RET"],
    "CPC": ["NKX2-5", "GATA6", "TBX5", "HAND2"],
    "Cardiac": ["TNNT2", "MYH6", "MYL7", "MYL2", "RYR2", "CACNA1C"],
}

# 1. ANNOTATED UMAPS, DOTPLOTS, & BALANCED WILCOXON DGE
for comp_name, h5_path in PAIRWISE_FILES.items():
    if not os.path.exists(h5_path):
        continue
    adata = sc.read_h5ad(h5_path)
    if "X_umap_15" in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_umap_15"].copy()

    id_col = "Manuscript_Identity" if "Manuscript_Identity" in adata.obs else "leiden_pairwise"

    # Global Annotated UMAP (Figs 3B, 5K, 6E)
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color=id_col, palette=cfg.ANNOTATION_COLORS, frameon=False, ax=ax, show=False)
    fig.savefig(os.path.join(cfg.FIG_DIR, f"{comp_name}_Global_Annotated_UMAP.pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig)

    # Lineage Annotation DotPlot (Figs S2A, S4A, S5B)
    pool = adata.raw.var_names if adata.raw is not None else adata.var_names
    seen = set()
    valid_markers = {}
    for k, genes in LINEAGE_MARKERS.items():
        g_list = [g for g in genes if g in pool and not (g in seen or seen.add(g))]
        if g_list:
            valid_markers[k] = g_list

    dp = sc.pl.DotPlot(adata, var_names=valid_markers, groupby=id_col, standard_scale="var", use_raw=(adata.raw is not None))
    dp.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.5)
    dp.savefig(os.path.join(cfg.FIG_DIR, f"{comp_name}_Lineage_Annotation_DotPlot.pdf"), bbox_inches="tight")

    # Balanced Wilcoxon DGE (Max 2,000 cells/cluster -> Table S1)
    np.random.seed(42)
    bal_idx = []
    for cl in adata.obs[id_col].unique():
        c_idx = adata.obs[adata.obs[id_col] == cl].index
        if len(c_idx) > 2000:
            c_idx = np.random.choice(c_idx, 2000, replace=False)
        bal_idx.extend(c_idx)
    adata_dge = adata[bal_idx].copy()
    sc.tl.rank_genes_groups(adata_dge, groupby=id_col, method="wilcoxon", use_raw=(adata_dge.raw is not None), pts=True)

    res = adata_dge.uns["rank_genes_groups"]
    dfs = []
    for g in res["names"].dtype.names:
        df_g = pd.DataFrame({
            "Cluster": g,
            "Gene": res["names"][g],
            "Log2FC": res["logfoldchanges"][g],
            "Adj_P_Value": res["pvals_adj"][g],
        })
        dfs.append(df_g[(df_g["Adj_P_Value"] < 0.05) & (df_g["Log2FC"] > 0.5)])
    pd.concat(dfs, ignore_index=True).to_csv(os.path.join(cfg.TABLE_DIR, f"{comp_name}_Wilcoxon_DGE_Markers.csv"), index=False)

# 2. HYPOXIA SPECIFIC ANALYSES (Figs 5M, 5N, 5O, S4C, S4E)
hypo_path = PAIRWISE_FILES["3_Hypoxia_vs_hSCA"]
if os.path.exists(hypo_path):
    adata_hyp = sc.read_h5ad(hypo_path)
    expr_hyp = adata_hyp.raw.to_adata() if adata_hyp.raw is not None else adata_hyp
    if expr_hyp.X.max() > 100:
        sc.pp.normalize_total(expr_hyp, target_sum=1e4)
        sc.pp.log1p(expr_hyp)

    expr_hyp.obs["Cluster_Condition"] = (
        expr_hyp.obs["Manuscript_Identity"].astype(str) + " (" +
        expr_hyp.obs["Biological_Sample"].astype(str).apply(lambda s: "Hypoxia" if "Hypo" in s else "Normoxia") + ")"
    )

    # Fig 5M: Apoptosis MatrixPlot
    apop_genes = [g for g in ["BAX", "CASP3", "CASP9", "BBC3"] if g in expr_hyp.var_names]
    target_conds = ["N1 (Normoxia)", "N1 (Hypoxia)", "CM1 (Normoxia)", "CM1 (Hypoxia)", "CM2 (Normoxia)", "CM2 (Hypoxia)"]
    sub_apop = expr_hyp[expr_hyp.obs["Cluster_Condition"].isin(target_conds)].copy()
    sub_apop.obs["Cluster_Condition"] = pd.Categorical(sub_apop.obs["Cluster_Condition"], categories=target_conds, ordered=True)

    mp_apop = sc.pl.matrixplot(sub_apop, apop_genes, groupby="Cluster_Condition", standard_scale="var", figsize=(4, 3.5), return_fig=True)
    mp_apop.style(edge_color="black", edge_lw=0.5)
    mp_apop.savefig(os.path.join(cfg.FIG_DIR, "Fig5M_Hypoxia_Apoptosis_Heatmap.pdf"), bbox_inches="tight")

    # Fig S4E: Hypoxia CM1 & CM2 Character Heatmap
    cm_conds = ["CM1 (Normoxia)", "CM1 (Hypoxia)", "CM2 (Normoxia)", "CM2 (Hypoxia)"]
    sub_cm_hyp = expr_hyp[expr_hyp.obs["Cluster_Condition"].isin(cm_conds)].copy()
    sub_cm_hyp.obs["Cluster_Condition"] = pd.Categorical(sub_cm_hyp.obs["Cluster_Condition"], categories=cm_conds, ordered=True)
    valid_phys = {k: [g for g in v if g in sub_cm_hyp.var_names] for k, v in cfg.PHYSIOLOGY_MODULES.items() if any(g in sub_cm_hyp.var_names for g in v)}
    mp_char = sc.pl.matrixplot(sub_cm_hyp, valid_phys, groupby="Cluster_Condition", standard_scale="var", cmap="viridis", swap_axes=True, return_fig=True)
    mp_char.style(edge_color="black", edge_lw=0.5)
    mp_char.savefig(os.path.join(cfg.FIG_DIR, "FigS4E_Hypoxia_CM_Character_Heatmap.pdf"), bbox_inches="tight")

    # Figs 5O & S4C: Per-Cluster Enrichr GO Biological Process 2023
    for cl in ["CM1", "CM2", "N1"]:
        sub_cl = expr_hyp[expr_hyp.obs["Manuscript_Identity"] == cl].copy()
        sc.tl.rank_genes_groups(sub_cl, groupby="Biological_Sample", groups=["hSCA06 Hypo"], reference="hSCA06 Normal", method="wilcoxon", use_raw=False)
        r_df = sc.get.rank_genes_groups_df(sub_cl, group="hSCA06 Hypo")
        up_genes = r_df[(r_df["pvals_adj"] < 0.05) & (r_df["logfoldchanges"] > 0.5)]["names"].head(150).tolist()
        if len(up_genes) >= 5:
            enr = gp.enrichr(gene_list=up_genes, gene_sets="GO_Biological_Process_2023", organism="human", outdir=None)
            enr.results.to_csv(os.path.join(cfg.TABLE_DIR, f"Fig5O_S4C_{cl}_Hypoxia_GO_BP.csv"), index=False)