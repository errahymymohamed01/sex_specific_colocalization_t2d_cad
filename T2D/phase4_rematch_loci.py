#!/usr/bin/env python3
"""
step3a_rematch_loci.py — Match ALL SNPs in ±500kb locus windows to eQTL data
=============================================================================

Purpose:
    eCAVIAR requires z-scores for ALL SNPs in a locus, not just the lead SNP.
    This script re-runs SNP matching from Step 2, but for every SNP within
    each ±500kb window around every lead variant.

Inputs:
    - GWAS lead SNPs:     GWAS/clumped/T2D_{sex}_lead_snps.tsv
    - GWAS summary stats: GWAS/T2D_{sex}_clean.tsv.gz  (comma-separated!)
    - GTEx eQTL data:     GTEx/{tissue}_clean.tsv       (tab-separated)

Outputs:
    Per-locus matched files:
        matched_loci/{sex}/{tissue}/locus_{lead_rsid}.tsv

    eCAVIAR z-score input files:
        ecaviar_input/{sex}/{tissue}/locus_{lead_rsid}.gwas.zscore
        ecaviar_input/{sex}/{tissue}/locus_{lead_rsid}.{gene_id}.eqtl.zscore

Author: Bioinformatics I Project — T2D GWAS Colocalization
Date:   2026-06-02
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import time
from pathlib import Path

# ── CONFIGURATION ────────────────────────────────────────────────────────
BASE_DIR = Path(r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\T2D\Cleaned1")

GWAS_DIR      = BASE_DIR / "GWAS"
CLUMPED_DIR   = GWAS_DIR / "clumped"
EQTL_DIR      = BASE_DIR / "GTEx"
MATCHED_DIR   = BASE_DIR / "matched_loci"
ECAVIAR_DIR   = BASE_DIR / "ecaviar_input"

SEXES = ["male"]

TISSUES = [
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Liver",
    "Pancreas",
]

WINDOW_KB = 500  # ±500kb around lead SNP

# Minimum SNPs in a locus to be useful for eCAVIAR
MIN_SNPS_PER_LOCUS = 2

# ── LOGGING ──────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(BASE_DIR / "scripts" / "step3a_rematch.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)


# ── ALLELE ALIGNMENT (same logic as snp_match.py) ────────────────────────
def check_allele_flip(gwas_ref, gwas_alt, eqtl_ref, eqtl_alt):
    """
    Check allele alignment between GWAS and eQTL.

    Returns:
        'match'    — alleles match directly
        'flip'     — alleles are swapped (need to negate eQTL z-score)
        'strand'   — A/T or C/G SNP (ambiguous strand — exclude)
        'mismatch' — alleles don't match at all (exclude)
    """
    ambiguous = {frozenset(["A", "T"]), frozenset(["C", "G"])}

    gwas_pair = frozenset([str(gwas_ref).upper(), str(gwas_alt).upper()])
    eqtl_pair = frozenset([str(eqtl_ref).upper(), str(eqtl_alt).upper()])

    # Check ambiguous strand
    if gwas_pair in ambiguous:
        return "strand"

    # Direct match
    if str(gwas_ref).upper() == str(eqtl_ref).upper() and str(gwas_alt).upper() == str(eqtl_alt).upper():
        return "match"

    # Flipped (effect allele swapped)
    if str(gwas_ref).upper() == str(eqtl_alt).upper() and str(gwas_alt).upper() == str(eqtl_ref).upper():
        return "flip"

    return "mismatch"


def load_gwas_summary(sex, loci_df):
    """Load full GWAS summary statistics (comma-separated .tsv.gz) in chunks and filter."""
    path = GWAS_DIR / f"T2D_{sex}_clean.tsv.gz"
    log.info(f"Loading GWAS summary stats: {path.name}")
    
    # Pre-compute valid intervals per chromosome for fast filtering
    valid_intervals = {}
    for _, row in loci_df.iterrows():
        chr_id = row['CHR']
        start_bp = row['BP'] - (WINDOW_KB * 1000)
        end_bp = row['BP'] + (WINDOW_KB * 1000)
        if chr_id not in valid_intervals:
            valid_intervals[chr_id] = []
        valid_intervals[chr_id].append((start_bp, end_bp))
        
    cols = ["rsid", "CHR", "BP", "reference_allele", "NEA", "beta", "se", "z_score", "p_value"]
    filtered_chunks = []
    chunksize = 2000000
    
    for chunk in pd.read_csv(path, sep=",", compression="gzip", usecols=cols, chunksize=chunksize):
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["BP"] = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "BP"])
        chunk["CHR"] = chunk["CHR"].astype(int)
        chunk["BP"] = chunk["BP"].astype(int)
        
        # Keep only rows where CHR and BP fall in valid_intervals
        mask = np.zeros(len(chunk), dtype=bool)
        for chr_id, intervals in valid_intervals.items():
            chr_mask = chunk['CHR'].values == chr_id
            for start, end in intervals:
                mask |= chr_mask & (chunk['BP'].values >= start) & (chunk['BP'].values <= end)
        filtered_chunks.append(chunk[mask])
        
    df = pd.concat(filtered_chunks, ignore_index=True)
    log.info(f"  Loaded {len(df):,} SNPs matching loci windows")

    return df


def load_lead_snps(sex):
    """Load lead SNPs from clumping results."""
    path = CLUMPED_DIR / f"T2D_{sex}_lead_snps.tsv"
    log.info(f"Loading lead SNPs: {path.name}")
    df = pd.read_csv(path, sep="\t")
    log.info(f"  {len(df)} lead SNPs loaded")

    # Ensure correct dtypes
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP"])
    df["CHR"] = df["CHR"].astype(int)
    df["BP"] = df["BP"].astype(int)

    return df


def load_eqtl(tissue):
    """Load eQTL data for a tissue."""
    path = EQTL_DIR / f"{tissue}_clean.tsv"
    log.info(f"  Loading eQTL: {path.name}")
    df = pd.read_csv(path, sep="\t")
    log.info(f"    {len(df):,} eQTL records loaded")

    # Ensure correct dtypes
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP"])
    df["CHR"] = df["CHR"].astype(int)
    df["BP"] = df["BP"].astype(int)

    # Compute eQTL z-score
    df["gtex_beta"] = pd.to_numeric(df["gtex_beta"], errors="coerce")
    df["gtex_se"] = pd.to_numeric(df["gtex_se"], errors="coerce")
    df["eqtl_zscore"] = df["gtex_beta"] / df["gtex_se"]

    return df


def get_locus_gwas_snps(gwas_df, lead_chr, lead_bp, window_kb=500):
    """Extract all GWAS SNPs within ±window_kb of the lead SNP."""
    window = window_kb * 1000
    mask = (
        (gwas_df["CHR"] == lead_chr)
        & (gwas_df["BP"] >= lead_bp - window)
        & (gwas_df["BP"] <= lead_bp + window)
    )
    return gwas_df[mask].copy()


def match_locus_snps(locus_gwas, eqtl_df, lead_rsid):
    """
    Match GWAS locus SNPs to eQTL by CHR+BP, perform allele alignment.

    Returns:
        matched_df: DataFrame with all matched SNPs and alignment info
    """
    # Merge on CHR + BP
    merged = pd.merge(
        locus_gwas,
        eqtl_df,
        on=["CHR", "BP"],
        how="inner",
        suffixes=("_gwas", "_eqtl"),
    )

    if len(merged) == 0:
        return pd.DataFrame()

    results = []
    for _, row in merged.iterrows():
        status = check_allele_flip(
            gwas_ref=row["reference_allele"],
            gwas_alt=row["NEA"],
            eqtl_ref=row["REF"],
            eqtl_alt=row["ALT"],
        )

        if status in ("mismatch", "strand"):
            continue

        # Flip eQTL z-score if alleles are swapped
        eqtl_z = row["eqtl_zscore"]
        eqtl_beta = row["gtex_beta"]
        if status == "flip":
            eqtl_z = -eqtl_z
            eqtl_beta = -eqtl_beta

        results.append(
            {
                "rsid": row["rsid"],
                "CHR": row["CHR"],
                "BP": row["BP"],
                "effect_allele": row["reference_allele"],
                "other_allele": row["NEA"],
                "gwas_beta": row["beta"],
                "gwas_se": row["se"],
                "gwas_zscore": row["z_score"],
                "gwas_pval": row["p_value"],
                "eqtl_beta": eqtl_beta,
                "eqtl_se": row["gtex_se"],
                "eqtl_zscore": eqtl_z,
                "eqtl_pval": row["gtex_pval"],
                "gene_id": row["gene_id"],
                "allele_status": status,
                "lead_rsid": lead_rsid,
            }
        )

    return pd.DataFrame(results)


def write_zscore_files(matched_df, lead_rsid, sex, tissue):
    """
    Write eCAVIAR z-score files for a locus.

    GWAS: one file with unique SNPs (deduplicated)
        locus_{lead_rsid}.gwas.zscore

    eQTL: one file per gene
        locus_{lead_rsid}.{gene_id}.eqtl.zscore
    """
    out_dir = ECAVIAR_DIR / sex / tissue
    out_dir.mkdir(parents=True, exist_ok=True)

    if len(matched_df) == 0:
        return 0, 0

    # ── GWAS z-score file: unique rsIDs only ─────────────────────────
    # If a SNP maps to multiple genes, keep only one GWAS z-score entry
    gwas_unique = matched_df.drop_duplicates(subset=["rsid"], keep="first")
    gwas_unique = gwas_unique.sort_values("BP")

    gwas_path = out_dir / f"locus_{lead_rsid}.gwas.zscore"
    gwas_unique[["rsid", "gwas_zscore"]].to_csv(
        gwas_path, sep=" ", index=False, header=False
    )

    # ── eQTL z-score files: one per gene ─────────────────────────────
    gene_count = 0
    for gene_id, gene_df in matched_df.groupby("gene_id"):
        # Deduplicate by rsid within each gene
        gene_unique = gene_df.drop_duplicates(subset=["rsid"], keep="first")
        gene_unique = gene_unique.sort_values("BP")

        # Only write if SNP count matches GWAS file
        # (eCAVIAR requires same SNP set in both files)
        # Use only SNPs present in GWAS z-score file
        gene_snps = gene_unique[gene_unique["rsid"].isin(gwas_unique["rsid"])]

        if len(gene_snps) < MIN_SNPS_PER_LOCUS:
            continue

        # Clean gene_id for filename (remove version suffix dots)
        gene_clean = gene_id.replace(".", "_")
        eqtl_path = out_dir / f"locus_{lead_rsid}.{gene_clean}.eqtl.zscore"

        # IMPORTANT: write in same SNP order as GWAS file
        gene_snps_ordered = gwas_unique[["rsid"]].merge(
            gene_snps[["rsid", "eqtl_zscore"]], on="rsid", how="inner"
        )

        gene_snps_ordered[["rsid", "eqtl_zscore"]].to_csv(
            eqtl_path, sep=" ", index=False, header=False
        )
        gene_count += 1

    return len(gwas_unique), gene_count


# ── MAIN ─────────────────────────────────────────────────────────────────
def main():
    log.info("=" * 65)
    log.info("Step 3a: Re-match ALL locus SNPs to eQTL data for eCAVIAR")
    log.info("=" * 65)

    start_time = time.time()
    overall_summary = []

    for sex in SEXES:
        log.info(f"\n{'='*60}")
        log.info(f"Processing sex: {sex}")
        log.info(f"{'='*60}")

        # Load lead SNPs first so we can filter
        lead_snps = load_lead_snps(sex)
        log.info(f"  {len(lead_snps)} lead SNPs to process")

        # Load full GWAS summary stats (filtered by loci to save memory)
        gwas_df = load_gwas_summary(sex, lead_snps)

        for tissue in TISSUES:
            log.info(f"\n  {'─'*50}")
            log.info(f"  Tissue: {tissue}")

            # Load eQTL data
            eqtl_path = EQTL_DIR / f"{tissue}_clean.tsv"
            if not eqtl_path.exists():
                log.warning(f"    eQTL file not found: {eqtl_path}")
                continue

            eqtl_df = load_eqtl(tissue)

            # Create output directories
            match_dir = MATCHED_DIR / sex / tissue
            match_dir.mkdir(parents=True, exist_ok=True)

            tissue_total_snps = 0
            tissue_total_loci = 0
            tissue_total_genes = 0
            tissue_skipped = 0

            for idx, lead in lead_snps.iterrows():
                lead_rsid = lead["rsid"]
                lead_chr = int(lead["CHR"])
                lead_bp = int(lead["BP"])

                # Get all GWAS SNPs in the ±500kb window
                locus_gwas = get_locus_gwas_snps(
                    gwas_df, lead_chr, lead_bp, window_kb=WINDOW_KB
                )

                if len(locus_gwas) == 0:
                    tissue_skipped += 1
                    continue

                # Match to eQTL
                matched = match_locus_snps(locus_gwas, eqtl_df, lead_rsid)

                if len(matched) < MIN_SNPS_PER_LOCUS:
                    tissue_skipped += 1
                    continue

                # Save per-locus matched file
                match_path = match_dir / f"locus_{lead_rsid}.tsv"
                matched.to_csv(match_path, sep="\t", index=False)

                # Write eCAVIAR z-score files
                n_snps, n_genes = write_zscore_files(matched, lead_rsid, sex, tissue)

                if n_snps >= MIN_SNPS_PER_LOCUS and n_genes > 0:
                    tissue_total_loci += 1
                    tissue_total_snps += n_snps
                    tissue_total_genes += n_genes

            log.info(
                f"    {tissue}: {tissue_total_loci} loci with >={MIN_SNPS_PER_LOCUS} "
                f"matched SNPs, {tissue_total_snps} total unique SNPs, "
                f"{tissue_total_genes} gene z-score files, "
                f"{tissue_skipped} skipped"
            )

            overall_summary.append(
                {
                    "sex": sex,
                    "tissue": tissue,
                    "loci_with_data": tissue_total_loci,
                    "total_matched_snps": tissue_total_snps,
                    "gene_zscore_files": tissue_total_genes,
                    "loci_skipped": tissue_skipped,
                }
            )

    # ── SUMMARY ──────────────────────────────────────────────────────
    elapsed = time.time() - start_time
    log.info(f"\n{'='*65}")
    log.info("SUMMARY")
    log.info(f"{'='*65}")

    summary_df = pd.DataFrame(overall_summary)
    log.info(f"\n{summary_df.to_string(index=False)}")

    # Save summary
    summary_path = BASE_DIR / "scripts" / "step3a_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    log.info(f"\nSummary saved to: {summary_path}")
    log.info(f"Elapsed time: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Report output directories
    log.info(f"\nOutput directories:")
    log.info(f"  Matched loci:    {MATCHED_DIR}")
    log.info(f"  eCAVIAR input:   {ECAVIAR_DIR}")

    return summary_df


if __name__ == "__main__":
    main()
