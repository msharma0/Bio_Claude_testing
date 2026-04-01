#!/usr/bin/env python3
"""
Improved Figure 6: Hallmarks of Aging expression heatmap.
Clearer labels, species grouping, better color scale, and annotation.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

RESULTS_DIR = "/Users/ambuj.mehta/code/Git_Repos/Bio_Claude_testing/results"
hallmark_csv = os.path.join(RESULTS_DIR, "Table3_Hallmarks_Aging_Expression.csv")

df = pd.read_csv(hallmark_csv)

# Convert Log2FC to numeric
df["Log2FC"] = pd.to_numeric(df["Log2FC"], errors="coerce")

# Dataset metadata for clearer labels
dataset_labels = {
    "GSE21942": "Multiple Sclerosis\n(Human, PBMCs)",
    "GSE8762":  "Huntington's Disease\n(Human, Lymphocytes)",
    "GSE35642": "Parkinsonism in vitro\n(Human, SK-N-MC cells)",
    "GSE23182": "Chronic Neurodegeneration\n(Mouse, Hippocampus)",
    "GSE26600": "MAM Genotoxin\n(Mouse, Brain)",
    "GSE47516": "Cstb KO / Epilepsy\n(Mouse, Cerebellum)",
    "GSE52118": "ALS Motor Pools\n(Mouse, Motor Neurons)",
    "GSE5786":  "PGC-1α KO / Huntington's\n(Mouse, Striatum)",
    "GSE63012": "IκBα KO / Neuroinflammation\n(Mouse, Hippocampus)",
    "GSE31458": "MPTP / Parkinson's\n(Mouse, PFC/Putamen)",
}

# Prettier hallmark names
hallmark_labels = {
    "Genomic_Instability": "Genomic\nInstability",
    "Telomere_Attrition": "Telomere\nAttrition",
    "Epigenetic_Alterations": "Epigenetic\nAlterations",
    "Loss_of_Proteostasis": "Loss of\nProteostasis",
    "Nutrient_Sensing": "Deregulated\nNutrient Sensing",
    "Mitochondrial_Dysfunction": "Mitochondrial\nDysfunction",
    "Cellular_Senescence": "Cellular\nSenescence",
    "Stem_Cell_Exhaustion": "Stem Cell\nExhaustion",
    "Altered_Intercellular_Comm": "Altered Intercell.\nCommunication",
}

# Compute mean Log2FC per hallmark per dataset
pivot = df.pivot_table(
    values="Log2FC",
    index="Dataset",
    columns="Hallmark",
    aggfunc="mean",
)

# Reorder: Human datasets first, then Mouse
human_order = ["GSE21942", "GSE8762", "GSE35642"]
mouse_order = ["GSE23182", "GSE26600", "GSE47516", "GSE52118", "GSE5786", "GSE63012", "GSE31458"]
row_order = [d for d in human_order + mouse_order if d in pivot.index]

# Reorder columns (hallmarks) in biological order
hallmark_order = [
    "Genomic_Instability", "Telomere_Attrition", "Epigenetic_Alterations",
    "Loss_of_Proteostasis", "Nutrient_Sensing", "Mitochondrial_Dysfunction",
    "Cellular_Senescence", "Stem_Cell_Exhaustion", "Altered_Intercellular_Comm",
]
col_order = [c for c in hallmark_order if c in pivot.columns]

pivot = pivot.loc[row_order, col_order]

# Rename for display
row_labels = [dataset_labels.get(d, d) for d in pivot.index]
col_labels = [hallmark_labels.get(c, c) for c in pivot.columns]

# Create figure
fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[20, 1], height_ratios=[3, 7],
                       hspace=0.05, wspace=0.15)

# Top panel: Human datasets
ax_human = fig.add_subplot(gs[0, 0])
human_data = pivot.loc[[d for d in human_order if d in pivot.index]]
human_row_labels = [dataset_labels.get(d, d) for d in human_data.index]

# Clamp extreme values for better color scaling
vmax = 1.5
vmin = -1.5

sns.heatmap(
    human_data.values.astype(float),
    xticklabels=[],  # will share with bottom
    yticklabels=human_row_labels,
    cmap="RdBu_r",
    center=0,
    vmin=vmin,
    vmax=vmax,
    annot=True,
    fmt=".2f",
    linewidths=1,
    linecolor="white",
    ax=ax_human,
    cbar=False,
    annot_kws={"fontsize": 10},
)
ax_human.set_title(
    "Hallmarks of Aging — Mean Gene Expression Change (Log2 Fold Change)\nin Neurodegenerative Disease vs Control",
    fontsize=14, fontweight="bold", pad=15,
)
ax_human.set_ylabel("HUMAN", fontsize=12, fontweight="bold", color="#c0392b")
ax_human.tick_params(axis="y", labelsize=10)

# Separator label
# Bottom panel: Mouse datasets
ax_mouse = fig.add_subplot(gs[1, 0])
mouse_data = pivot.loc[[d for d in mouse_order if d in pivot.index]]
mouse_row_labels = [dataset_labels.get(d, d) for d in mouse_data.index]

sns.heatmap(
    mouse_data.values.astype(float),
    xticklabels=col_labels,
    yticklabels=mouse_row_labels,
    cmap="RdBu_r",
    center=0,
    vmin=vmin,
    vmax=vmax,
    annot=True,
    fmt=".2f",
    linewidths=1,
    linecolor="white",
    ax=ax_mouse,
    cbar=False,
    annot_kws={"fontsize": 10},
)
ax_mouse.set_ylabel("MOUSE", fontsize=12, fontweight="bold", color="#2980b9")
ax_mouse.tick_params(axis="y", labelsize=10)
ax_mouse.set_xticklabels(ax_mouse.get_xticklabels(), rotation=45, ha="right", fontsize=10)

# Colorbar
ax_cbar = fig.add_subplot(gs[:, 1])
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = matplotlib.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=ax_cbar)
cbar.set_label("Mean Log2 Fold Change\n(Disease vs Control)", fontsize=11)
cbar.ax.tick_params(labelsize=10)

# Add explanatory text at the bottom
fig.text(
    0.5, 0.01,
    "Each cell = average Log2FC of representative genes for that hallmark.  "
    "Red = upregulated in disease/KO.  Blue = downregulated.  "
    "Values clamped to ±1.5 for display (actual extremes noted in Table 3 CSV).",
    ha="center", fontsize=9, fontstyle="italic", color="gray",
)

plt.savefig(os.path.join(RESULTS_DIR, "Fig6_Hallmarks_Aging_Heatmap.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(RESULTS_DIR, "Fig6_Hallmarks_Aging_Heatmap.pdf"), bbox_inches="tight")
plt.close()
print("Saved improved Fig6!")
