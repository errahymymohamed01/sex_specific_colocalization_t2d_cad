import requests
import pandas as pd
import numpy as np
import time
import os
from io import StringIO
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================
TOKEN      = "13b50f3ce4aa"
LDLINK     = "https://ldlink.nih.gov/LDlinkRest"
TRAIT      = "CAD"
MATCHED_DIR= r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\matched"
OUTPUT_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\ld_matrices"
POP        = "CEU"
WINDOW_KB  = 50    # 100kb window → ~100 SNPs → matrice ~100×100
BATCH      = 20     # max SNPs par appel LDmatrix
MAX_RETRY  = 8
# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Collecter tous les lead SNPs uniques depuis les fichiers matched
all_leads = set()
for fname in os.listdir(MATCHED_DIR):
    if fname.endswith("_matched.tsv"):
        df = pd.read_csv(os.path.join(MATCHED_DIR, fname), sep="\t")
        all_leads.update(df["rsid"].tolist())

lead_snps = sorted(all_leads)
print(f"Lead SNPs to process: {len(lead_snps)}")
print(f"Window: ±{WINDOW_KB}kb | Batch: {BATCH} | Expected matrix size: ~100×100")
print(f"Estimated time: ~{len(lead_snps) * 2} minutes\n")


def ldproxy(snp, pop, token, window_kb):
    """Get all SNPs in ±window_kb around lead SNP."""
    try:
        r = requests.get(f"{LDLINK}/ldproxy", params={
            "var":    snp,
            "pop":    pop,
            "r2_d":   "r2",
            "window": window_kb * 1000,
            "token":  token,
        }, timeout=120)

        df  = pd.read_csv(StringIO(r.text), sep="\t")
        col = df.columns[0]
        snps = df[df[col].astype(str).str.startswith("rs")][col].tolist()
        result = list(dict.fromkeys([snp] + snps))
        return result
    except Exception as e:
        print(f"    LDproxy error: {e.__class__.__name__}")
        return [snp]


def ldmatrix_batch(snp_list, pop, token):
    """Compute LD matrix for a batch of SNPs."""
    for attempt in range(MAX_RETRY):
        try:
            r = requests.post(
                f"{LDLINK}/ldmatrix?token={token}",
                data={"snps": "\n".join(snp_list), "pop": pop, "r2_d": "r2"},
                timeout=180
            )
            if r.status_code == 200 and not r.text.startswith("{"):
                return r.text
            print(f"    HTTP {r.status_code} — {r.text[:60].strip()}")
        except Exception as e:
            print(f"    Connection error: {e.__class__.__name__}")
        wait = 5 * (attempt + 1)
        print(f"    Retry {attempt+1}/{MAX_RETRY} — wait {wait}s")
        time.sleep(wait)
    return None


def compute_ld_matrix(lead_snp, pop, token, window_kb, batch_size):
    """
    Full pipeline for one locus:
    1. LDproxy → get all SNPs in window
    2. LDmatrix in batches → assemble full matrix
    """
    all_snps = ldproxy(lead_snp, pop, token, window_kb)
    print(f"    LDproxy: {len(all_snps)} SNPs in ±{window_kb}kb window")

    if len(all_snps) < 2:
        return None, []

    n      = len(all_snps)
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, 1.0)  # diagonal = 1

    chunks_done = 0
    chunks_total = (n - 1) // batch_size + 1

    for i in range(0, n, batch_size):
        chunk = all_snps[i:i + batch_size]
        chunks_done += 1
        print(f"    Chunk {chunks_done}/{chunks_total} ({len(chunk)} SNPs)...")

        raw = ldmatrix_batch(chunk, pop, token)
        if raw is None:
            print(f"    Chunk {chunks_done} failed — skipping")
            continue

        try:
            sub_df = pd.read_csv(StringIO(raw), sep="\t", index_col=0)
            sub_df = sub_df.apply(pd.to_numeric, errors="coerce").fillna(0)

            for ri, rsnp in enumerate(sub_df.index):
                if rsnp in all_snps:
                    gi = all_snps.index(rsnp)
                    for ci, csnp in enumerate(sub_df.columns):
                        if csnp in all_snps:
                            gj = all_snps.index(csnp)
                            matrix[gi, gj] = sub_df.iloc[ri, ci]
        except Exception as e:
            print(f"    Parse error: {e.__class__.__name__}")

        time.sleep(5)  # pause entre chunks

    return matrix, all_snps


def save_ld_matrix(matrix, snps, out_prefix):
    ld_path  = Path(f"{out_prefix}.ld")
    snp_path = Path(f"{out_prefix}.snplist")
    np.savetxt(ld_path, matrix, fmt="%.6f", delimiter="\t")
    snp_path.write_text("\n".join(snps))
    return ld_path


# ============================================================
# MAIN PIPELINE
# ============================================================
ok = fail = skipped = 0

for i, lead in enumerate(lead_snps, 1):
    out_prefix = os.path.join(OUTPUT_DIR, lead)
    ld_path    = Path(f"{out_prefix}.ld")

    # Resume: skip if already computed and not empty
    if ld_path.exists() and ld_path.stat().st_size > 100:
        size = ld_path.stat().st_size / 1024
        print(f"[{i:3d}/{len(lead_snps)}] SKIP {lead} (already done — {size:.0f} KB)")
        skipped += 1
        ok += 1
        continue

    print(f"\n[{i:3d}/{len(lead_snps)}] Processing {lead}...")

    try:
        matrix, snps = compute_ld_matrix(
            lead_snp  = lead,
            pop       = POP,
            token     = TOKEN,
            window_kb = WINDOW_KB,
            batch_size= BATCH
        )
    except Exception as e:
        print(f"    CRASH — {e.__class__.__name__}: {e}")
        fail += 1
        time.sleep(5)
        continue

    if matrix is not None and len(snps) > 1:
        save_ld_matrix(matrix, snps, out_prefix)
        size = ld_path.stat().st_size / 1024
        print(f"    OK — {len(snps)}×{len(snps)} matrix — {size:.0f} KB")
        ok += 1
    else:
        print(f"    FAIL — could not compute matrix for {lead}")
        fail += 1

    time.sleep(5)  # pause entre loci

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*50}")
print(f"RÉSUMÉ FINAL — {TRAIT} LD Matrices")
print(f"{'='*50}")
print(f"  OK      : {ok}/{len(lead_snps)}")
print(f"  Skipped : {skipped}")
print(f"  Failed  : {fail}")
print(f"\nMatrices saved in: {OUTPUT_DIR}")

# Show matrix sizes
print(f"\nMatrix sizes:")
for lead in lead_snps:
    ld_path = Path(os.path.join(OUTPUT_DIR, f"{lead}.ld"))
    snp_path = Path(os.path.join(OUTPUT_DIR, f"{lead}.snplist"))
    if ld_path.exists() and snp_path.exists():
        n_snps = len(snp_path.read_text().strip().split("\n"))
        size   = ld_path.stat().st_size / 1024
        print(f"  {lead}: {n_snps}×{n_snps} — {size:.0f} KB")
    else:
        print(f"  {lead}: MISSING")
