#!/usr/bin/env python3
"""
Step 6+7: Compare female vs male CLPP + generate all plots
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

# ── Load results ──────────────────────────────────────────────────────────────
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--results", required=True, help="clpp_all_results.csv")
parser.add_argument("--outdir",  default="ecaviar_plots")
args = parser.parse_args()

Path(args.outdir).mkdir(parents=True, exist_ok=True)
df = pd.read_csv(args.results)
print(f"Loaded {len(df)} combinations")

THRESHOLD = 0.01
STRONG    = 0.10

# ── Colors ────────────────────────────────────────────────────────────────────
COL_F = "#E07B8A"   # pink for female
COL_M = "#5B8DB8"   # blue for male
COL_S = "#6A9E78"   # green for shared

TISSUES = [
    "Adipose_Subcutaneous","Adipose_Visceral_Omentum",
    "Artery_Aorta","Artery_Coronary","Artery_Tibial",
    "Liver","Pancreas"
]
TISSUE_SHORT = {
    "Adipose_Subcutaneous":       "Adipose\nSubcut.",
    "Adipose_Visceral_Omentum":   "Adipose\nVisceral",
    "Artery_Aorta":               "Artery\nAorta",
    "Artery_Coronary":            "Artery\nCoronary",
    "Artery_Tibial":              "Artery\nTibial",
    "Liver":                      "Liver",
    "Pancreas":                   "Pancreas",
}

female_df = df[df["sex"]=="female"]
male_df   = df[df["sex"]=="male"]

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1: CLPP distribution by sex (violin)
# ═══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7,5))
parts = ax.violinplot(
    [female_df["clpp"].values, male_df["clpp"].values],
    positions=[1,2], showmedians=True, showextrema=True
)
for pc, col in zip(parts["bodies"], [COL_F, COL_M]):
    pc.set_facecolor(col); pc.set_alpha(0.7)
parts["cmedians"].set_color("black")
ax.set_xticks([1,2]); ax.set_xticklabels(["Female","Male"], fontsize=13)
ax.set_ylabel("CLPP score", fontsize=12)
ax.set_title("CLPP Distribution by Sex", fontsize=14, fontweight="bold")
ax.axhline(THRESHOLD, color="red", linestyle="--", alpha=0.6, label="CLPP=0.01")
ax.axhline(STRONG,    color="orange", linestyle="--", alpha=0.6, label="CLPP=0.10")
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f"{args.outdir}/1_clpp_violin.png", dpi=150)
plt.close()
print("Plot 1: violin done")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2: Colocalized loci count by tissue and sex (grouped bar)
# ═══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10,5))
x = np.arange(len(TISSUES))
w = 0.35
f_counts = [female_df[(female_df["tissue"]==t)&(female_df["clpp"]>=THRESHOLD)].shape[0] for t in TISSUES]
m_counts = [  male_df[(  male_df["tissue"]==t)&(  male_df["clpp"]>=THRESHOLD)].shape[0] for t in TISSUES]
ax.bar(x - w/2, f_counts, w, label="Female", color=COL_F, alpha=0.85)
ax.bar(x + w/2, m_counts, w, label="Male",   color=COL_M, alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels([TISSUE_SHORT[t] for t in TISSUES], fontsize=10)
ax.set_ylabel("# Colocalized loci (CLPP ≥ 0.01)", fontsize=11)
ax.set_title("Colocalized Loci by Tissue and Sex", fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(f"{args.outdir}/2_coloc_by_tissue.png", dpi=150)
plt.close()
print("Plot 2: bar chart done")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3: Sex-specific vs shared colocalization (stacked bar per tissue)
# ═══════════════════════════════════════════════════════════════════════════════
# For loci present in BOTH sexes at the same locus+tissue+gene
shared_key = ["locus","tissue","gene"]
f_sig = female_df[female_df["clpp"]>=THRESHOLD].set_index(shared_key)["clpp"]
m_sig =   male_df[  male_df["clpp"]>=THRESHOLD].set_index(shared_key)["clpp"]

female_only = len(f_sig.index.difference(m_sig.index))
male_only   = len(m_sig.index.difference(f_sig.index))
shared      = len(f_sig.index.intersection(m_sig.index))

fig, axes = plt.subplots(1,2, figsize=(12,5))

# Pie chart
axes[0].pie(
    [female_only, male_only, shared],
    labels=[f"Female only\n({female_only})", f"Male only\n({male_only})", f"Shared\n({shared})"],
    colors=[COL_F, COL_M, COL_S],
    autopct="%1.1f%%", startangle=140,
    textprops={"fontsize":11}
)
axes[0].set_title("Colocalization Specificity\n(CLPP ≥ 0.01)", fontsize=13, fontweight="bold")

# Per tissue stacked bar
tissue_stats = []
for t in TISSUES:
    ft = set(female_df[(female_df["tissue"]==t)&(female_df["clpp"]>=THRESHOLD)]
             .apply(lambda r: (r["locus"],r["gene"]), axis=1))
    mt = set(  male_df[(  male_df["tissue"]==t)&(  male_df["clpp"]>=THRESHOLD)]
             .apply(lambda r: (r["locus"],r["gene"]), axis=1))
    tissue_stats.append({
        "tissue": TISSUE_SHORT[t],
        "female_only": len(ft - mt),
        "shared":      len(ft & mt),
        "male_only":   len(mt - ft),
    })
ts = pd.DataFrame(tissue_stats)
x = np.arange(len(ts))
axes[1].bar(x, ts["female_only"], label="Female only", color=COL_F, alpha=0.85)
axes[1].bar(x, ts["shared"],      bottom=ts["female_only"], label="Shared", color=COL_S, alpha=0.85)
axes[1].bar(x, ts["male_only"],   bottom=ts["female_only"]+ts["shared"], label="Male only", color=COL_M, alpha=0.85)
axes[1].set_xticks(x); axes[1].set_xticklabels(ts["tissue"], fontsize=9)
axes[1].set_ylabel("# Colocalized locus-gene pairs", fontsize=11)
axes[1].set_title("Sex Specificity by Tissue", fontsize=13, fontweight="bold")
axes[1].legend(fontsize=10)
plt.tight_layout()
plt.savefig(f"{args.outdir}/3_sex_specificity.png", dpi=150)
plt.close()
print("Plot 3: specificity done")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4: Heatmap — top loci × tissues (CLPP, female vs male side by side)
# ═══════════════════════════════════════════════════════════════════════════════
top_loci = (df[df["clpp"]>=STRONG]
            .groupby("locus")["clpp"].max()
            .sort_values(ascending=False)
            .head(25).index.tolist())

def make_heatmap_matrix(sub_df, loci, tissues):
    mat = pd.DataFrame(index=loci, columns=tissues, dtype=float).fillna(0)
    for _, r in sub_df.iterrows():
        if r["locus"] in loci and r["tissue"] in tissues:
            mat.loc[r["locus"], r["tissue"]] = max(mat.loc[r["locus"], r["tissue"]], r["clpp"])
    return mat

f_mat = make_heatmap_matrix(female_df, top_loci, TISSUES)
m_mat = make_heatmap_matrix(  male_df, top_loci, TISSUES)

fig, axes = plt.subplots(1, 2, figsize=(14, 9), sharey=True)
kw = dict(vmin=0, vmax=0.5, cmap="YlOrRd", linewidths=0.3, linecolor="grey")

sns.heatmap(f_mat, ax=axes[0], **kw, cbar=False,
            xticklabels=[TISSUE_SHORT[t].replace("\n"," ") for t in TISSUES])
axes[0].set_title("Female", fontsize=14, fontweight="bold", color=COL_F)
axes[0].set_xlabel(""); axes[0].set_ylabel("Lead SNP", fontsize=11)

im = sns.heatmap(m_mat, ax=axes[1], **kw,
                 xticklabels=[TISSUE_SHORT[t].replace("\n"," ") for t in TISSUES])
axes[1].set_title("Male", fontsize=14, fontweight="bold", color=COL_M)
axes[1].set_ylabel("")

for ax in axes:
    ax.tick_params(axis="x", rotation=35, labelsize=9)
    ax.tick_params(axis="y", rotation=0,  labelsize=8)

fig.suptitle("CLPP Heatmap — Top 25 Loci × Tissues", fontsize=15, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(f"{args.outdir}/4_clpp_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Plot 4: heatmap done")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 5: Female vs Male CLPP scatter (shared locus-tissue-gene combos)
# ═══════════════════════════════════════════════════════════════════════════════
f_idx = female_df.set_index(shared_key)["clpp"]
m_idx =   male_df.set_index(shared_key)["clpp"]
both  = f_idx.index.intersection(m_idx.index)

if len(both) >= 5:
    fc = f_idx.loc[both].values
    mc = m_idx.loc[both].values
    fig, ax = plt.subplots(figsize=(7,6))
    sc = ax.scatter(fc, mc, alpha=0.5, c=fc+mc, cmap="viridis", s=30, edgecolors="none")
    plt.colorbar(sc, ax=ax, label="Sum CLPP")
    ax.axline((0,0), slope=1, color="grey", linestyle="--", alpha=0.5, label="y=x")
    ax.axhline(THRESHOLD, color=COL_M, linestyle=":", alpha=0.6)
    ax.axvline(THRESHOLD, color=COL_F, linestyle=":", alpha=0.6)
    ax.set_xlabel("Female CLPP", fontsize=12, color=COL_F)
    ax.set_ylabel("Male CLPP",   fontsize=12, color=COL_M)
    ax.set_title(f"Female vs Male CLPP\n({len(both)} shared locus-tissue-gene pairs)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    corr = np.corrcoef(fc, mc)[0,1]
    ax.text(0.05, 0.92, f"r = {corr:.3f}", transform=ax.transAxes, fontsize=11,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"))
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/5_female_vs_male_scatter.png", dpi=150)
    plt.close()
    print(f"Plot 5: scatter done (r={corr:.3f})")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary table
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(f"Total combinations tested:         {len(df)}")
print(f"Colocalized CLPP>=0.01:            {(df['clpp']>=0.01).sum()}")
print(f"Strongly colocalized CLPP>=0.10:   {(df['clpp']>=0.10).sum()}")
print(f"\nFemale: {(female_df['clpp']>=0.01).sum()} colocalized")
print(f"Male:   {(  male_df['clpp']>=0.01).sum()} colocalized")
print(f"\nSex-specific breakdown (CLPP>=0.01):")
print(f"  Female only: {female_only}")
print(f"  Male only:   {male_only}")
print(f"  Shared:      {shared}")
print(f"\nTop 10 hits:")
print(df[df["clpp"]>=STRONG].nlargest(10,"clpp")[["sex","locus","tissue","gene","clpp"]].to_string(index=False))
print(f"\nPlots saved to: {args.outdir}/")
