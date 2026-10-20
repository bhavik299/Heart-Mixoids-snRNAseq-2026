"""
02_Pairwise_Reintegration_and_Clustree.py
=========================================
1. Partitions the global master atlas into targeted pairwise comparisons.
2. Wipes global reductions and recomputes local batch-aware HVGs, PCA, and Harmony.
3. Evaluates Leiden resolution stability (0.2 - 1.2) via Clustree and Euclidean Silhouette scores.
4. Enforces the strict upper multiplet cutoff (< 10,000 detected genes).
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

master_path = os.path.join(cfg.QC_DIR, "hSCA_Canonical_Integrated_Master.h5ad")
master_adata = sc.read_h5ad(master_path)

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
OPT_RES = 0.2

def compute_silhouette_curve(adata, resolutions, rep="X_pca_harmony"):
    X = adata.obsm[rep]
    rows = []
    for res in resolutions:
        labels = adata.obs[f"leiden_{res}"].values
        n_clust = len(np.unique(labels))
        if n_clust < 2:
            rows.append({"resolution": res, "n_clusters": n_clust, "silhouette_mean": np.nan})
            continue
        idx = np.random.default_rng(0).choice(adata.n_obs, min(adata.n_obs, 5000), replace=False)
        s = silhouette_score(X[idx], labels[idx], metric="euclidean", sample_size=min(len(idx), 3000), random_state=0)
        rows.append({"resolution": res, "n_clusters": n_clust, "silhouette_mean": s})
    return pd.DataFrame(rows)

for comp_name, samples in cfg.PAIRWISE_COMPARISONS.items():
    sub = master_adata[master_adata.obs["Biological_Sample"].isin(samples)].copy()
    if sub.raw is not None:
        sub = sub.raw.to_adata()
    elif "counts" in sub.layers:
        sub.X = sub.layers["counts"].copy()

    for k in list(sub.obsm.keys()): del sub.obsm[k]
    for k in list(sub.varm.keys()): del sub.varm[k]
    for k in list(sub.obsp.keys()): del sub.obsp[k]

    # Enforce Upper QC Limit (< 10,000 genes)
    sc.pp.calculate_qc_metrics(sub, percent_top=None, log1p=False, inplace=True)
    sub = sub[sub.obs["n_genes_by_counts"] < cfg.MAX_GENES_UPPER].copy()

    # Local Normalization & Harmony Reintegration
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.raw = sub

    sc.pp.highly_variable_genes(sub, n_top_genes=3000, flavor="seurat", batch_key="Sample")
    sub_hvg = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(sub_hvg, max_value=10)
    sc.tl.pca(sub_hvg, svd_solver="arpack", n_comps=40)
    sub.obsm["X_pca"] = sub_hvg.obsm["X_pca"].copy()

    hm_out = hm.run_harmony(sub.obsm["X_pca"], sub.obs, ["Sample"])
    z_corr = np.array(hm_out.Z_corr)
    sub.obsm["X_pca_harmony"] = z_corr.T if z_corr.shape[0] != sub.n_obs else z_corr

    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=40)
    sc.tl.umap(sub)
    sub.obsm["X_umap_15"] = sub.obsm["X_umap"].copy()

    for res in RESOLUTIONS:
        sc.tl.leiden(sub, resolution=res, key_added=f"leiden_{res}")
    sub.obs["leiden_pairwise"] = sub.obs[f"leiden_{OPT_RES}"].copy()

    # Save Silhouette Diagnostics
    sil_df = compute_silhouette_curve(sub, RESOLUTIONS)
    sil_df.to_csv(os.path.join(cfg.TABLE_DIR, f"{comp_name}_Silhouette_Scores.csv"), index=False)

    out_h5ad = os.path.join(cfg.PAIRWISE_DIR, f"{comp_name[2:]}_Raw1.h5ad")
    sub.write(out_h5ad)
    print(f"Saved pairwise manifold: {out_h5ad}")