import pandas as pd
import os

# ── CONFIG ────────────────────────────────────────────────
GWAS_DIR  = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\clumped"
EQTL_DIR  = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\GTEx"
OUT_DIR   = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\matched"
os.makedirs(OUT_DIR, exist_ok=True)

TRAIT = "CAD"
SEXES = ["female", "male"]

EQTL_FILES = [
    "Adipose_Subcutaneous_clean.tsv",
    "Adipose_Visceral_Omentum_clean.tsv",
    "Artery_Aorta_clean.tsv",
    "Artery_Coronary_clean.tsv",
    "Artery_Tibial_clean.tsv",
    "Liver_clean.tsv",
    "Pancreas_clean.tsv"
]
# ─────────────────────────────────────────────────────────


def check_allele_flip(gwas_ref, gwas_alt, eqtl_ref, eqtl_alt):
    ambiguous = {frozenset(['A', 'T']), frozenset(['C', 'G'])}
    gwas_ref  = str(gwas_ref).upper().strip()
    gwas_alt  = str(gwas_alt).upper().strip()
    eqtl_ref  = str(eqtl_ref).upper().strip()
    eqtl_alt  = str(eqtl_alt).upper().strip()
    if frozenset([gwas_ref, gwas_alt]) in ambiguous:
        return 'strand'
    if gwas_ref == eqtl_ref and gwas_alt == eqtl_alt:
        return 'match'
    if gwas_ref == eqtl_alt and gwas_alt == eqtl_ref:
        return 'flip'
    return 'mismatch'


def match_gwas_eqtl(gwas_df, eqtl_df, sex, tissue):
    eqtl_df = eqtl_df.copy()
    eqtl_df['eqtl_zscore'] = eqtl_df['gtex_beta'] / eqtl_df['gtex_se']

    # Harmonize types
    gwas_df = gwas_df.copy()
    gwas_df['CHR'] = gwas_df['CHR'].astype(str)
    gwas_df['BP']  = gwas_df['BP'].astype(int)
    eqtl_df['CHR'] = eqtl_df['CHR'].astype(str)
    eqtl_df['BP']  = eqtl_df['BP'].astype(int)

    merged = pd.merge(gwas_df, eqtl_df, on=['CHR', 'BP'], how='inner')
    print(f"    SNPs matched by position: {len(merged):,}")
    if len(merged) == 0:
        return pd.DataFrame()

    results = []
    for _, row in merged.iterrows():
        if str(row['NEA']).strip() == '.':
            continue

        status = check_allele_flip(
            gwas_ref=row['reference_allele'],
            gwas_alt=row['NEA'],
            eqtl_ref=row['REF'],
            eqtl_alt=row['ALT']
        )
        if status in ('mismatch', 'strand'):
            continue

        eqtl_z = -row['eqtl_zscore'] if status == 'flip' else row['eqtl_zscore']

        results.append({
            'rsid'          : row['rsid'],
            'CHR'           : row['CHR'],
            'BP'            : row['BP'],
            'effect_allele' : row['reference_allele'],
            'other_allele'  : row['NEA'],
            'gwas_beta'     : row['beta'],
            'gwas_se'       : row['se'],
            'gwas_zscore'   : row['z_score'],
            'gwas_pval'     : row['p_value'],
            'eqtl_beta'     : row['gtex_beta'],
            'eqtl_se'       : row['gtex_se'],
            'eqtl_zscore'   : eqtl_z,
            'eqtl_pval'     : row['gtex_pval'],
            'gene_id'       : row['gene_id'],
            'allele_status' : status,
            'sex'           : sex,
            'tissue'        : tissue
        })

    return pd.DataFrame(results)


# ── MAIN ──────────────────────────────────────────────────
summary = []

for sex in SEXES:
    print(f"\n{'='*50}")
    print(f"Sex: {sex}")

    # ← SAME AS T2D: use lead_snps.tsv (one SNP per locus)
    gwas_path = os.path.join(GWAS_DIR, f"{TRAIT}_{sex}_lead_snps.tsv")
    if not os.path.exists(gwas_path):
        print(f"  WARNING: not found -> {gwas_path}")
        continue

    gwas = pd.read_csv(gwas_path, sep='\t')
    gwas['CHR'] = gwas['CHR'].astype(str)
    gwas['BP']  = gwas['BP'].astype(int)
    print(f"  GWAS lead SNPs: {len(gwas):,}")

    for eqtl_file in EQTL_FILES:
        tissue    = eqtl_file.replace("_clean.tsv", "")
        eqtl_path = os.path.join(EQTL_DIR, eqtl_file)
        print(f"\n  Tissue: {tissue}")

        if not os.path.exists(eqtl_path):
            print(f"    WARNING: not found -> {eqtl_path}")
            continue

        eqtl = pd.read_csv(eqtl_path, sep='\t', low_memory=False)
        print(f"    eQTL SNPs loaded: {len(eqtl):,}")

        matched = match_gwas_eqtl(gwas, eqtl, sex, tissue)

        if len(matched) == 0:
            print(f"    No matches found")
            continue

        n_match = (matched['allele_status'] == 'match').sum()
        n_flip  = (matched['allele_status'] == 'flip').sum()
        print(f"    After allele check: {len(matched):,} SNPs")
        print(f"    Direct match: {n_match} | Flipped: {n_flip}")

        out_file = os.path.join(OUT_DIR, f"{TRAIT}_{sex}_{tissue}_matched.tsv")
        matched.to_csv(out_file, sep='\t', index=False)
        print(f"    Saved -> {TRAIT}_{sex}_{tissue}_matched.tsv")

        summary.append({
            'sex'         : sex,
            'tissue'      : tissue,
            'matched_snps': len(matched),
            'flipped'     : n_flip
        })

print(f"\n{'='*50}")
print("SUMMARY:")
summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))
summary_df.to_csv(os.path.join(OUT_DIR, "matching_summary.tsv"), sep='\t', index=False)
print("\nDone.")
