import pandas as pd
import os

GWAS_DIR   = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned"
OUTPUT_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\clumped"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TRAIT       = "CAD"
SEXES       = ["female", "male"]
PVAL_THRESH = 5e-8
CLUMP_KB    = 500

def load_gwas(trait, sex):
    fname = f"{trait}_{sex}_clean.tsv"
    fpath = os.path.join(GWAS_DIR, fname)
    if not os.path.exists(fpath):
        print(f"  WARNING: not found -> {fpath}")
        return None, None

    print(f"  Reading {fname} in chunks...")
    sig_chunks = []
    total = 0

    for chunk in pd.read_csv(fpath, sep=",", chunksize=200_000):
        total += len(chunk)
        sig_chunk = chunk[chunk["p_value"] < PVAL_THRESH]
        if len(sig_chunk) > 0:
            sig_chunks.append(sig_chunk)

    sig = pd.concat(sig_chunks, ignore_index=True) if sig_chunks else pd.DataFrame()
    print(f"  Total SNPs: {total:,} | Significant (p < 5e-8): {len(sig):,}")
    return sig, fpath   # return sig only + path for window step

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

def define_windows(lead_snps_df, gwas_fpath, window_kb=500):
    """Read the full GWAS file in chunks and collect SNPs within lead SNP windows."""
    # Build a set of (CHR, min_bp, max_bp) intervals for fast filtering
    intervals = [
        (row["CHR"], row["BP"] - window_kb * 1000, row["BP"] + window_kb * 1000, row["rsid"], row["BP"], row["p_value"])
        for _, row in lead_snps_df.iterrows()
    ]

    loci_chunks = []
    for chunk in pd.read_csv(gwas_fpath, sep=",", chunksize=200_000):
        matched = []
        for chrom, bp_min, bp_max, lead_rsid, lead_pos, lead_pval in intervals:
            win = chunk[
                (chunk["CHR"] == chrom) &
                (chunk["BP"]  >= bp_min) &
                (chunk["BP"]  <= bp_max)
            ].copy()
            if len(win) > 0:
                win["lead_rsid"] = lead_rsid
                win["lead_pos"]  = lead_pos
                win["lead_pval"] = lead_pval
                matched.append(win)
        if matched:
            loci_chunks.append(pd.concat(matched, ignore_index=True))

    if loci_chunks:
        return pd.concat(loci_chunks, ignore_index=True).drop_duplicates()
    return pd.DataFrame()

summary = []
for sex in SEXES:
    print(f"\n==================================================")
    print(f"Processing: {TRAIT} — {sex}")
    sig, gwas_fpath = load_gwas(TRAIT, sex)
    if sig is None or len(sig) == 0:
        print("  No significant SNPs — skipping")
        continue

    print("  Clumping...")
    lead_snps = distance_clump(sig, window_kb=CLUMP_KB)
    print(f"  Independent loci: {len(lead_snps):,}")

    print("  Defining locus windows (chunked read)...")
    loci_df = define_windows(lead_snps, gwas_fpath, window_kb=CLUMP_KB)
    print(f"  Total SNPs in windows: {len(loci_df):,}")

    lead_snps.to_csv(os.path.join(OUTPUT_DIR, f"{TRAIT}_{sex}_lead_snps.tsv"), sep="\t", index=False)
    loci_df.to_csv(os.path.join(OUTPUT_DIR, f"{TRAIT}_{sex}_loci.tsv"), sep="\t", index=False)
    print(f"  Saved -> {TRAIT}_{sex}_lead_snps.tsv")
    print(f"  Saved -> {TRAIT}_{sex}_loci.tsv")

    summary.append({"trait": TRAIT, "sex": sex, "sig_snps": len(sig), "loci": len(lead_snps)})

print(f"\n==================================================")
print("DONE — SUMMARY:")
summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))
summary_df.to_csv(os.path.join(OUTPUT_DIR, "clumping_summary.tsv"), sep="\t", index=False)