"""
05_YAP1_MEIS1_Continuum_and_Signaling.py
========================================
1. Assigns purified CMs into the 4 discrete YAP1/MEIS1 continuum states on raw counts:
   - State 0: YAP1- / MEIS1- (Double Negative)
   - State 1: YAP1+ / MEIS1- (Proliferative)
   - State 2: YAP1+ / MEIS1+ (Mature)
   - State 3: YAP1- / MEIS1+ (Mature-Primed)
2. Exports fractional penetrance CSVs for ERBB2, ERBB4, YAP1, MEIS1, and RYR2 (Figs 3N, 3O, 6H, 6N, 6O, S3G).
3. Generates 4-state and 3-step maturation heatmaps and ERBB2/4 dotplots (Figs 3P, 3Q, S3I, S5C, S5D).
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import config as cfg

warnings.filterwarnings("ignore")
cfg.apply_publication_aesthetics()

DATASETS = {
    "hSCM_vs_CM14t": os.path.join(cfg.PAIRWISE_DIR, "hSCA_vs_CM14t_Pristine_CM_Manifold.h5ad"),
    "hPCM_vs_CM14t": os.path.join(cfg.PAIRWISE_DIR, "hPSCA_vs_CM14t_Pure_Thresholded.h5ad"),
}

def get_expr_vector(adata_obj, gene):
    if gene not in adata_obj.var_names:
        return np.zeros(adata_obj.n_obs)
    arr = adata_obj[:, gene].X
    return arr.toarray().flatten() if hasattr(arr, "toarray") else arr.flatten()

for label, path in DATASETS.items():
    if not os.path.exists(path):
        continue
    adata_full = sc.read_h5ad(path)
    adata_raw = adata_full.raw.to_adata() if adata_full.raw is not None else adata_full.copy()

    yap1 = get_expr_vector(adata_raw, "YAP1")
    meis1 = get_expr_vector(adata_raw, "MEIS1")

    states = np.full(adata_raw.n_obs, "0. Double Negative (YAP1-/MEIS1-)", dtype=object)
    states[(yap1 > 0) & (meis1 <= 0)] = "1. Proliferative (YAP1+/MEIS1-)"
    states[(yap1 > 0) & (meis1 > 0)]  = "2. Mature (YAP1+/MEIS1+)"
    states[(yap1 <= 0) & (meis1 > 0)] = "3. Mature-Primed (YAP1-/MEIS1+)"

    state_order = [
        "0. Double Negative (YAP1-/MEIS1-)",
        "1. Proliferative (YAP1+/MEIS1-)",
        "2. Mature (YAP1+/MEIS1+)",
        "3. Mature-Primed (YAP1-/MEIS1+)",
    ]
    adata_raw.obs["Continuum_State"] = pd.Categorical(states, categories=state_order, ordered=True)

    # Export Sample-Resolved State Proportions (Figs 3O, 6O, S3G)
    prop_df = pd.crosstab(adata_raw.obs["Biological_Sample"], adata_raw.obs["Continuum_State"], normalize="index") * 100
    prop_df.to_csv(os.path.join(cfg.TABLE_DIR, f"{label}_YAP1_MEIS1_State_Percentages.csv"))

    # Export Fractional Penetrance of Key Signaling Markers (Figs 3N, 6H, 6N)
    pen_rows = []
    for samp in adata_raw.obs["Biological_Sample"].unique():
        sub_s = adata_raw[adata_raw.obs["Biological_Sample"] == samp]
        for g in ["ERBB2", "ERBB4", "YAP1", "MEIS1", "RYR2"]:
            if g in sub_s.var_names:
                pct = (get_expr_vector(sub_s, g) > 0).mean() * 100
                pen_rows.append({"Dataset": label, "Sample": samp, "Gene": g, "Percent_Positive": round(pct, 2)})
    pd.DataFrame(pen_rows).to_csv(os.path.join(cfg.TABLE_DIR, f"{label}_Key_Marker_Penetrance.csv"), index=False)

    # Render Sample-Specific 4-State Maturation Heatmaps (Figs 3P, S3I, S5D)
    valid_mods = {k: [g for g in v if g in adata_raw.var_names] for k, v in cfg.PHYSIOLOGY_MODULES.items() if any(g in adata_raw.var_names for g in v)}
    for samp in adata_raw.obs["Biological_Sample"].unique():
        sub_s = adata_raw[adata_raw.obs["Biological_Sample"] == samp].copy()
        safe_s = str(samp).replace(" ", "_")
        mp = sc.pl.matrixplot(sub_s, var_names=valid_mods, groupby="Continuum_State", standard_scale="var", swap_axes=True, cmap="viridis", return_fig=True)
        mp.style(edge_color="black", edge_lw=0.5)
        mp.savefig(os.path.join(cfg.FIG_DIR, f"{label}_{safe_s}_4State_Continuum_Heatmap.pdf"), bbox_inches="tight")

    # Render 3-Step ERBB2/ERBB4 Receptor DotPlot (Figs 3Q, S5C)
    sub_3step = adata_raw[adata_raw.obs["Continuum_State"] != "0. Double Negative (YAP1-/MEIS1-)"].copy()
    sub_3step.obs["Continuum_State"] = sub_3step.obs["Continuum_State"].cat.remove_unused_categories()
    erbb_genes = [g for g in ["ERBB2", "ERBB4"] if g in sub_3step.var_names]
    if erbb_genes:
        dp = sc.pl.dotplot(sub_3step, var_names=erbb_genes, groupby="Continuum_State", standard_scale="var", swap_axes=True, cmap="Blues", return_fig=True)
        dp.style(cmap="Blues", dot_edge_color="black", dot_edge_lw=0.5, grid=True)
        dp.savefig(os.path.join(cfg.FIG_DIR, f"{label}_3Step_ERBB2_4_DotPlot.pdf"), bbox_inches="tight")