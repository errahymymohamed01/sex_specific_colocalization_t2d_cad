import argparse, os
import pandas as pd
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--gwas_female",  required=True)
parser.add_argument("--gwas_male",    required=True)
parser.add_argument("--gtex_dir",     required=True)
parser.add_argument("--ld_dir",       required=True)
parser.add_argument("--lead_female",  required=True)
parser.add_argument("--lead_male",    required=True)
parser.add_argument("--outdir",       default="ecaviar_inputs")
parser.add_argument("--window",       type=int, default=500000)
parser.add_argument("--min_snps",     type=int, default=3)
args = parser.parse_args()

WINDOW = args.window

def to_int_chr(val):
    try:
        return int(val)
    except:
        return None

# Load lead SNPs
leads = {
    "female": pd.read_csv(args.lead_female, sep="\t").dropna(subset=["rsid","CHR","BP"]),
    "male":   pd.read_csv(args.lead_male,   sep="\t").dropna(subset=["rsid","CHR","BP"]),
}

# Build windows per sex: dict of {chrom: [(pos_min, pos_max, rsid), ...]}
def build_windows(lead_df, window):
    wins = {}
    for _, r in lead_df.iterrows():
        c = int(r["CHR"]); p = int(r["BP"])
        wins.setdefault(c, []).append((p - window, p + window, str(r["rsid"])))
    return wins

windows = {
    "female": build_windows(leads["female"], WINDOW),
    "male":   build_windows(leads["male"],   WINDOW),
}

# Load GWAS in chunks - only keep SNPs inside any locus window
print("Loading GWAS (chunked to save memory)...")
gwas = {}
for sex, fpath in [("female", args.gwas_female), ("male", args.gwas_male)]:
    wins = windows[sex]
    kept = []
    chunk_n = 0
    for chunk in pd.read_csv(fpath, sep=",", low_memory=False, chunksize=100000,
                              usecols=["rsid","CHR","BP","z_score"]):
        chunk_n += 1
        chunk = chunk.dropna(subset=["z_score"])
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR"])
        chunk["CHR"] = chunk["CHR"].astype(int)
        chunk["BP"]  = chunk["BP"].astype(int)
        # Keep only rows inside a locus window
        mask = pd.Series(False, index=chunk.index)
        for chrom, win_list in wins.items():
            chr_mask = chunk["CHR"] == chrom
            if not chr_mask.any():
                continue
            for lo, hi, _ in win_list:
                mask |= (chr_mask & (chunk["BP"] >= lo) & (chunk["BP"] <= hi))
        kept.append(chunk[mask])
        if chunk_n % 10 == 0:
            print(f"  {sex}: {chunk_n*100000:,} rows scanned, {sum(len(k) for k in kept):,} kept", end="\r")
    gwas[sex] = pd.concat(kept, ignore_index=True)
    print(f"  {sex}: done — {len(gwas[sex]):,} SNPs in windows          ")

# Load GTEx tissue files
print("\nLoading GTEx tissue files...")
gtex = {}
for f in Path(args.gtex_dir).glob("*_clean.tsv"):
    tissue = f.stem.replace("_clean", "")
    df = pd.read_csv(f, sep="\t", low_memory=False,
                     usecols=["CHR","BP","gene_id","gtex_beta","gtex_se"])
    df["CHR"] = df["CHR"].apply(to_int_chr)
    df = df.dropna(subset=["CHR"])
    df["CHR"] = df["CHR"].astype(int)
    df["BP"]  = df["BP"].astype(int)
    df["eqtl_zscore"] = df["gtex_beta"] / df["gtex_se"]
    gtex[tissue] = df
    print(f"  {tissue}: {len(df):,} rows")

# Process loci
summary = []
for sex in ["female", "male"]:
    lead_df = leads[sex]
    gwas_df = gwas[sex]
    ld_sex  = Path(args.ld_dir) / sex
    print(f"\n{'='*50}\nProcessing {sex.upper()} loci...\n{'='*50}")

    for _, locus in lead_df.iterrows():
        rsid  = str(locus["rsid"])
        chrom = int(locus["CHR"])
        pos   = int(locus["BP"])

        snps_file = ld_sex / f"{rsid}.snps"
        ld_file   = ld_sex / f"{rsid}.ld"
        if not snps_file.exists() or not ld_file.exists():
            continue

        ld_snps = open(snps_file).read().strip().split("\n")
        ld_mat  = np.loadtxt(ld_file)
        if ld_mat.ndim == 0: ld_mat = np.array([[float(ld_mat)]])
        if ld_mat.ndim == 1: ld_mat = ld_mat.reshape(1,1)
        snp_idx = {s: i for i, s in enumerate(ld_snps)}

        gwas_win = gwas_df[
            (gwas_df["CHR"] == chrom) &
            (gwas_df["BP"] >= pos - WINDOW) &
            (gwas_df["BP"] <= pos + WINDOW)
        ].copy()

        for tissue, gtex_df in gtex.items():
            gtex_win = gtex_df[
                (gtex_df["CHR"] == chrom) &
                (gtex_df["BP"] >= pos - WINDOW) &
                (gtex_df["BP"] <= pos + WINDOW)
            ].copy()
            if gtex_win.empty: continue

            merged = gwas_win.merge(gtex_win[["CHR","BP","gene_id","eqtl_zscore"]],
                                    on=["CHR","BP"], how="inner")
            if merged.empty: continue
            merged = merged[merged["rsid"].isin(snp_idx)].copy()
            if len(merged) < args.min_snps: continue

            for gene in merged["gene_id"].unique():
                gene_df = merged[merged["gene_id"] == gene].copy()
                if len(gene_df) < args.min_snps: continue

                gene_df["ld_order"] = gene_df["rsid"].map(snp_idx)
                gene_df = gene_df.sort_values("ld_order").drop_duplicates("rsid")
                kept_snps = gene_df["rsid"].tolist()
                kept_idx  = [snp_idx[s] for s in kept_snps]
                sub_ld    = ld_mat[np.ix_(kept_idx, kept_idx)]

                gene_short = gene.split(".")[0]
                out_dir = Path(args.outdir) / sex / tissue / gene_short
                out_dir.mkdir(parents=True, exist_ok=True)
                prefix = out_dir / rsid

                gene_df[["rsid","z_score"]].to_csv(f"{prefix}_gwas.z", sep="\t", index=False, header=False)
                gene_df[["rsid","eqtl_zscore"]].to_csv(f"{prefix}_eqtl.z", sep="\t", index=False, header=False)
                np.savetxt(f"{prefix}.ld", sub_ld, fmt="%.6f", delimiter=" ")
                summary.append({"sex":sex,"locus":rsid,"tissue":tissue,"gene":gene_short,"n_snps":len(kept_snps)})

        print(f"  {rsid}: {len([s for s in summary if s['locus']==rsid])} tissue-gene combos")

print(f"\n{'='*50}")
print(f"Total combinations: {len(summary)}")
summary_df = pd.DataFrame(summary)
if len(summary_df) > 0:
    print(summary_df.groupby(["sex","tissue"])["locus"].count().to_string())
    summary_df.to_csv(Path(args.outdir) / "ecaviar_input_summary.csv", index=False)
    print("Summary saved to ecaviar_input_summary.csv")
