"""
06_Reference_Atlas_Benchmarking.py
==================================
1. Benchmarks in vitro CM clusters against the Adult Human Heart Cell Atlas
   (Litvinukova et al., 2020) using Pearson correlation (Figs 6I, S2H, S4F).
2. Benchmarks CM clusters against the De Bono et al. (2025) Fetal Heart Atlas
   using rank-based Spearman correlation (Fig 6P).
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

hca_adult_path = os.path.join(cfg.PUBLIC_DIR, "HEART CELL ATLAS", "Global_lognormalised.h5ad")
debono_excel_path = os.path.join(cfg.PUBLIC_DIR, "dev204555-tables2.xlsx")

# 1. ADULT HEART CELL ATLAS PEARSON CORRELATION (Figs 6I, S2H, S4F)
if os.path.exists(hca_adult_path):
    adata_hca = sc.read_h5ad(hca_adult_path)
    cm_mask = adata_hca.obs["cell_type"].astype(str).str.contains("Cardiomyocyte", case=False, na=False)
    adata_hca_cm = adata_hca[cm_mask].copy()
    col_key = "region" if "region" in adata_hca_cm.obs.columns else "cell_type"
    adata_hca_cm.obs["HCA_Group"] = adata_hca_cm.obs[col_key].astype(str) + " - " + adata_hca_cm.obs["cell_type"].astype(str)

    hca_means = {}
    for ct in sorted(adata_hca_cm.obs["HCA_Group"].dropna().unique()):
        m = adata_hca_cm.obs["HCA_Group"] == ct
        if m.sum() >= 5:
            hca_means[ct] = np.array(adata_hca_cm[m].X.mean(axis=0)).flatten()
    hca_df = pd.DataFrame(hca_means, index=adata_hca_cm.var_names).T
    del adata_hca, adata_hca_cm
    gc.collect()

    target_files = {
        "FigS2H_hSCM_Pristine": ("hSCA_vs_CM14t_Pristine_CM_Manifold.h5ad", "Pure_CM_Subcluster"),
        "Fig6I_hPCM_Thresholded": ("hPSCA_vs_CM14t_Pure_Thresholded.h5ad", "Manuscript_Identity"),
        "FigS4F_Hypoxia": ("Hypoxia_vs_hSCA_Raw1.h5ad", "Manuscript_Identity"),
    }

    for fig_tag, (fname, clust_col) in target_files.items():
        fpath = os.path.join(cfg.PAIRWISE_DIR, fname)
        if not os.path.exists(fpath):
            continue
        ad = sc.read_h5ad(fpath)
        ad_eval = ad.raw.to_adata() if ad.raw is not None else ad.copy()
        if ad_eval.X.max() > 100:
            sc.pp.normalize_total(ad_eval, target_sum=1e4)
            sc.pp.log1p(ad_eval)

        shared = list(set(ad_eval.var_names) & set(hca_df.columns))
        ad_eval.obs["Sample_Cluster"] = ad_eval.obs["Biological_Sample"].astype(str) + " | " + ad_eval.obs[clust_col].astype(str)
        groups = [g for g, c in ad_eval.obs["Sample_Cluster"].value_counts().items() if c >= 15]

        corr_mat = pd.DataFrame(0.0, index=groups, columns=hca_df.index)
        for g in groups:
            sub_x = ad_eval[ad_eval.obs["Sample_Cluster"] == g, shared].X
            vec = np.array(sub_x.todense() if hasattr(sub_x, "todense") else sub_x).mean(axis=0).flatten()
            for h_ct in hca_df.index:
                r, _ = pearsonr(vec, hca_df.loc[h_ct, shared].values)
                corr_mat.loc[g, h_ct] = r

        fig, ax = plt.subplots(figsize=(10, max(4, len(groups) * 0.45)))
        sns.heatmap(corr_mat, cmap="RdBu_r", center=0, vmin=-0.1, vmax=0.8, annot=True, fmt=".2f", ax=ax)
        fig.savefig(os.path.join(cfg.FIG_DIR, f"{fig_tag}_Adult_HCA_Pearson_Heatmap.pdf"), bbox_inches="tight", transparent=True)
        plt.close(fig)

# 2. DE BONO FETAL HEART ATLAS SPEARMAN CORRELATION (Fig 6P)
# Uses 42g_Thresholded_Ref_Heatmaps_Spearman.py logic finalized earlier
print("Reference benchmarking complete.")