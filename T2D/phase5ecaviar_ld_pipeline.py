#!/usr/bin/env python3
"""
eCAVIAR LD Matrix Pipeline using LDlink API
Computes LD matrices for T2D female and male lead SNPs

Usage:
    python ecaviar_ld_pipeline.py --token YOUR_LDLINK_TOKEN

Get your free LDlink token at: https://ldlink.nih.gov/?tab=apiaccess
"""

import argparse
import os
import time
import requests
import pandas as pd
import numpy as np
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────────────────
WINDOW_KB    = 500        # ±500 kb window around each lead SNP
POPULATION   = "CEU+TSI+FIN+GBR+IBS"   # 1000G EUR super-population
R2_THRESHOLD = 0.0        # keep all SNPs (eCAVIAR needs full matrix)
LDLINK_URL   = "https://ldlink.nih.gov/LDlinkRest"
SLEEP_SEC    = 1.0        # pause between API calls (be polite to NIH server)

# ── Argument parsing ────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="eCAVIAR LD matrix pipeline via LDlink")
parser.add_argument("--token",   required=True,  help="2066f4ff77f9")
parser.add_argument("--female",  default="T2D_female_lead_snps.tsv")
parser.add_argument("--male",    default="T2D_male_lead_snps.tsv")
parser.add_argument("--outdir",  default="ecaviar_ld_matrices")
parser.add_argument("--window",  type=int, default=WINDOW_KB,
                    help="Window size in kb (default 500)")
parser.add_argument("--pop",     default=POPULATION,
                    help="1000G population code(s)")
args = parser.parse_args()

Path(args.outdir).mkdir(parents=True, exist_ok=True)
female_dir = Path(args.outdir) / "female"
male_dir   = Path(args.outdir) / "male"
female_dir.mkdir(exist_ok=True)
male_dir.mkdir(exist_ok=True)


# ── Helper functions ────────────────────────────────────────────────────────────

def get_ld_matrix(snp_list, population, token):
    """
    Call LDmatrix endpoint for a list of rsIDs.
    Returns a tuple: (snp_list_returned, r2_matrix_df) or (None, None) on error.
    """
    payload = {
        "snps": "\n".join(snp_list),
        "pop":  population,
        "r2_d": "r2",
        "token": token,
    }
    try:
        r = requests.post(f"{LDLINK_URL}/ldmatrix", data=payload, timeout=120)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        if len(lines) < 2 or lines[0].startswith("error"):
            print(f"    LDlink error: {lines[0][:120]}")
            return None, None
        # Parse TSV: first column = row label, first row = header
        from io import StringIO
        df = pd.read_csv(StringIO(r.text), sep="\t", index_col=0)
        return df.columns.tolist(), df
    except Exception as e:
        print(f"    Request failed: {e}")
        return None, None


def get_proxy_snps_in_window(lead_snp, chrom, pos_bp, window_kb, population, token):
    """
    Use LDproxy to get all SNPs in LD within the window around a lead SNP.
    Returns list of rsIDs (including the lead SNP).
    """
    payload = {
        "var":   lead_snp,
        "pop":   population,
        "r2_d":  "r2",
        "window": window_kb * 1000,
        "token": token,
    }
    try:
        r = requests.get(f"{LDLINK_URL}/ldproxy", params=payload, timeout=120)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        if lines[0].startswith("error") or len(lines) < 2:
            print(f"    LDproxy error for {lead_snp}: {lines[0][:100]}")
            return [lead_snp]   # fallback: just the lead SNP itself
        from io import StringIO
        df = pd.read_csv(StringIO(r.text), sep="\t")
        if "RS_Number" in df.columns:
            snps = df[df["RS_Number"].str.startswith("rs", na=False)]["RS_Number"].tolist()
            return list(dict.fromkeys([lead_snp] + snps))  # deduplicated, lead SNP first
        return [lead_snp]
    except Exception as e:
        print(f"    LDproxy request failed for {lead_snp}: {e}")
        return [lead_snp]


def save_ecaviar_ld(lead_snp, snp_list, r2_df, out_dir):
    """
    Save LD matrix in eCAVIAR format:
      - <lead_snp>.ld  : space-separated R² matrix (no header, no row names)
      - <lead_snp>.snps: one rsID per line (same order as matrix rows/cols)
    """
    matrix = r2_df.values
    np.fill_diagonal(matrix, 1.0)          # ensure diagonal = 1
    matrix = np.clip(matrix, 0, 1)         # clip any floating-point noise

    ld_path   = Path(out_dir) / f"{lead_snp}.ld"
    snps_path = Path(out_dir) / f"{lead_snp}.snps"

    np.savetxt(ld_path,   matrix, fmt="%.6f", delimiter=" ")
    with open(snps_path, "w") as fh:
        fh.write("\n".join(snp_list) + "\n")

    return ld_path, snps_path


def process_locus(row, out_dir, population, token, window_kb):
    lead_snp = row["rsid"]
    chrom    = int(row["CHR"])
    pos_bp   = int(row["BP"])

    ld_path = Path(out_dir) / f"{lead_snp}.ld"
    if ld_path.exists():
        print(f"  [SKIP] {lead_snp} already done")
        return True

    print(f"  Processing {lead_snp}  chr{chrom}:{pos_bp:,}")

    # Step 1: get SNPs in window
    snp_list = get_proxy_snps_in_window(lead_snp, chrom, pos_bp, window_kb, population, token)
    time.sleep(SLEEP_SEC)

    if len(snp_list) < 2:
        print(f"    Only {len(snp_list)} SNP(s) in window — writing 1×1 matrix")
        dummy = pd.DataFrame([[1.0]], index=[lead_snp], columns=[lead_snp])
        save_ecaviar_ld(lead_snp, [lead_snp], dummy, out_dir)
        return True

    print(f"    {len(snp_list)} SNPs in ±{window_kb}kb window")

    # Step 2: get LD matrix
    # LDmatrix has a limit of ~1000 SNPs per call; chunk if needed
    MAX_SNPS = 900
    if len(snp_list) > MAX_SNPS:
        print(f"    Truncating to {MAX_SNPS} SNPs (lead SNP kept first)")
        snp_list = snp_list[:MAX_SNPS]

    returned_snps, r2_df = get_ld_matrix(snp_list, population, token)
    time.sleep(SLEEP_SEC)

    if r2_df is None:
        print(f"    LD matrix failed for {lead_snp} — skipping")
        return False

    # Step 3: save
    ld_path, snps_path = save_ecaviar_ld(lead_snp, returned_snps, r2_df, out_dir)
    print(f"    Saved {len(returned_snps)}×{len(returned_snps)} matrix → {ld_path.name}")
    return True


# ── Main ────────────────────────────────────────────────────────────────────────

def run(snp_file, out_dir, label):
    df = pd.read_csv(snp_file, sep="\t")
    df = df.dropna(subset=["rsid", "CHR", "BP"])
    df["CHR"] = df["CHR"].astype(int)
    df["BP"]  = df["BP"].astype(int)

    print(f"\n{'='*60}")
    print(f"  {label}: {len(df)} loci")
    print(f"  Output: {out_dir}")
    print(f"{'='*60}")

    success = 0
    for _, row in df.iterrows():
        ok = process_locus(row, out_dir, args.pop, args.token, args.window)
        if ok:
            success += 1

    print(f"\n  Done: {success}/{len(df)} loci processed successfully")


if __name__ == "__main__":
    run(args.female, female_dir, "FEMALE")
    run(args.male,   male_dir,   "MALE")
    print("\nAll done! LD matrices saved to:", args.outdir)
    print("\nNext step — run eCAVIAR per locus:")
    print("  eCAVIAR -l <locus>.ld -z <locus>.z -o <locus>_ecaviar")
