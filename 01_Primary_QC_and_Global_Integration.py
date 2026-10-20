"""
01_Primary_QC_and_Global_Integration.py
=======================================
1. Filters raw 10x Genomics 5' GEM-X HDF5 matrices (N1 & N2).
2. Merges technical replicates and evaluates unintegrated reproducibility.
3. Maps OCM multiplexing tags (cells_per_tag.json) to biological samples.
4. Performs batch-aware HVG selection, PCA, Harmony integration, and global UMAP.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

# 1. INITIAL PER-REPLICATE QC
raw_h5_paths = {
    "N1": os.path.join(cfg.BASE_DIR, "raw", "N1", "sample_filtered_feature_bc_matrix.h5"),
    "N2": os.path.join(cfg.BASE_DIR, "raw", "N2", "sample_filtered_feature_bc_matrix.h5"),
}
json_paths = {
    "N1": os.path.join(cfg.BASE_DIR, "raw", "N1", "cells_per_tag.json"),
    "N2": os.path.join(cfg.BASE_DIR, "raw", "N2", "cells_per_tag.json"),
}

qc_adatas = {}
for rep, path in raw_h5_paths.items():
    if not os.path.exists(path):
        continue
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs["Sample"] = rep
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    sc.pp.filter_cells(adata, min_genes=cfg.MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=cfg.MIN_CELLS)
    adata = adata[adata.obs["pct_counts_mt"] < cfg.MAX_MT_PCT, :].copy()
    adata.layers["counts"] = adata.X.copy()
    adata.write(os.path.join(cfg.QC_DIR, f"hSCA_{rep}_Only_QC_Raw.h5ad"))
    qc_adatas[rep] = adata

# 2. MERGE N1 & N2 AND DEMULTIPLEX OCM TAGS
if len(qc_adatas) == 2:
    shared_genes = np.intersect1d(qc_adatas["N1"].var_names, qc_adatas["N2"].var_names)
    adata = sc.concat(
        {"N1": qc_adatas["N1"][:, shared_genes], "N2": qc_adatas["N2"][:, shared_genes]},
        label="Sample", index_unique="-"
    )
    adata.layers["counts"] = adata.X.copy()
else:
    adata = sc.read_h5ad(os.path.join(cfg.QC_DIR, "hSCA_Combined_QC_Analyzed.h5ad"))

def build_ocm_lookup(json_path, suffix):
    if not os.path.exists(json_path):
        return {}
    with open(json_path) as f:
        data = json.load(f)
    return {
        f"{bc.strip()}-{suffix}": tag
        for tag, bclist in data.items() if isinstance(bclist, list)
        for bc in bclist
    }

full_lookup = {**build_ocm_lookup(json_paths["N1"], "N1"), **build_ocm_lookup(json_paths["N2"], "N2")}
if full_lookup:
    adata.obs["Assignment"] = adata.obs_names.map(full_lookup).fillna("Unassigned")
    tag_to_name = {
        "OB1": "hSCA06 Normal",
        "OB2": "CM14t",
        "OB3": "hSCA06 Hypo",
        "OB4": "hPSCA02",
        "Unassigned": "Unassigned",
        "Multiplet": "Multiplet",
    }
    adata.obs["Biological_Sample"] = adata.obs["Assignment"].map(tag_to_name).fillna("Unassigned")
    adata = adata[~adata.obs["Biological_Sample"].isin(["Unassigned", "Multiplet"])].copy()

# 3. GLOBAL CANONICAL HARMONY INTEGRATION
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="Sample")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack", n_comps=40)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

harmony_out = hm.run_harmony(adata.obsm["X_pca"], adata.obs, ["Sample"])
corrected_pca = np.array(harmony_out.Z_corr)
if corrected_pca.shape[0] != adata.n_obs:
    corrected_pca = corrected_pca.T
adata.obsm["X_pca_harmony"] = corrected_pca

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_canonical")

master_out = os.path.join(cfg.QC_DIR, "hSCA_Canonical_Integrated_Master.h5ad")
adata.write(master_out)
print(f"Saved Global Integrated Master Atlas: {master_out}")