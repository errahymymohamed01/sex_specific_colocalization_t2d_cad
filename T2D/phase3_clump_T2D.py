import pandas as pd
import os

GWAS_DIR   = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\T2D\Cleaned\GWAS"
OUTPUT_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\T2D\Cleaned\GWAS\clumped"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TRAITS      = ["T2D"]
SEXES       = ["female", "male"]
PVAL_THRESH = 5e-8
CLUMP_KB    = 500

def load_gwas(trait, sex):
    fname = f"{trait}_{sex}_clean.tsv.gz"
    fpath = os.path.join(GWAS_DIR, fname)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, sep=",", compression="gzip")
        print(f"  Loaded {fname} — {len(df):,} SNPs")
        print(f"  Columns: {df.columns.tolist()}")
        return df
    print(f"  WARNING: not found -> {fpath}")
    return None

def distance_clump(sig_df, window_kb=500):
    sig = sig_df.sort_values("p_value").copy()
    kept = []
    excluded = set()
    for _, row in sig.iterrows():
        snp = row["rsid"]
        if snp in excluded:
            continue
        kept.append(row)
        same_chr = sig[sig["CHR"] == row["CHR"]]
        nearby = same_chr[
            (same_chr["BP"] >= row["BP"] - window_kb * 1000) &
            (same_chr["BP"] <= row["BP"] + window_kb * 1000)
        ]
        excluded.update(nearby["rsid"].tolist())
    return pd.DataFrame(kept).reset_index(drop=True)

def define_windows(lead_snps_df, full_gwas_df, window_kb=500):
    loci = []
    for _, lead in lead_snps_df.iterrows():
        window = full_gwas_df[
            (full_gwas_df["CHR"] == lead["CHR"]) &
            (full_gwas_df["BP"]  >= lead["BP"] - window_kb * 1000) &
            (full_gwas_df["BP"]  <= lead["BP"] + window_kb * 1000)
        ].copy()
        window["lead_rsid"] = lead["rsid"]
        window["lead_pos"]  = lead["BP"]
        window["lead_pval"] = lead["p_value"]
        loci.append(window)
    if loci:
        return pd.concat(loci, ignore_index=True)
    return pd.DataFrame()

summary = []
for trait in TRAITS:
    for sex in SEXES:
        print(f"\n==================================================")
        print(f"Processing: {trait} — {sex}")
        df = load_gwas(trait, sex)
        if df is None:
            continue
        sig = df[df["p_value"] < PVAL_THRESH].copy()
        print(f"  Significant SNPs (p < 5e-8): {len(sig):,}")
        if len(sig) == 0:
            print("  No significant SNPs — skipping")
            continue
        print("  Clumping...")
        lead_snps = distance_clump(sig, window_kb=CLUMP_KB)
        print(f"  Independent loci: {len(lead_snps):,}")
        loci_df = define_windows(lead_snps, df, window_kb=CLUMP_KB)
        print(f"  Total SNPs in windows: {len(loci_df):,}")
        lead_snps.to_csv(os.path.join(OUTPUT_DIR, f"{trait}_{sex}_lead_snps.tsv"), sep="\t", index=False)
        loci_df.to_csv(os.path.join(OUTPUT_DIR, f"{trait}_{sex}_loci.tsv"), sep="\t", index=False)
        print(f"  Saved -> {trait}_{sex}_lead_snps.tsv")
        print(f"  Saved -> {trait}_{sex}_loci.tsv")
        summary.append({"trait": trait, "sex": sex, "sig_snps": len(sig), "loci": len(lead_snps)})

print(f"\n==================================================")
print("DONE — SUMMARY:")
summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))
summary_df.to_csv(os.path.join(OUTPUT_DIR, "clumping_summary.tsv"), sep="\t", index=False)
