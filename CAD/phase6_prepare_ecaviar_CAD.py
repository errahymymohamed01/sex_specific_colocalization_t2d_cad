"""
Phase 6a — Prépare les fichiers .z et .ld pour eCAVIAR (CAD)
Crée la structure: ecaviar_inputs/sex/tissue/gene/locus_gwas.z + locus_eqtl.z + locus.ld
"""
import pandas as pd
import numpy as np
import os
from pathlib import Path

# ── CONFIG ────────────────────────────────────────────────
MATCHED_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\matched"
LD_DIR      = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\ld_matrices"
OUT_DIR     = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\ecaviar_inputs"
SUMMARY_OUT = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned\ecaviar_summary.csv"
TRAIT       = "CAD"
SEXES       = ["female", "male"]
TISSUES     = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
    "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
    "Liver", "Pancreas"
]
# ─────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)
summary_rows = []

for sex in SEXES:
    for tissue in TISSUES:
        mfile = os.path.join(MATCHED_DIR, f"{TRAIT}_{sex}_{tissue}_matched.tsv")
        if not os.path.exists(mfile):
            continue

        df = pd.read_csv(mfile, sep="\t")
        if df.empty:
            continue

        # Group by locus (rsid = lead SNP) and gene
        for (locus, gene), grp in df.groupby(["rsid", "gene_id"]):
            ld_path = os.path.join(LD_DIR, sex, f"{locus}.ld")
            snp_path = os.path.join(LD_DIR, sex, f"{locus}.snplist")

            # Output directory for this combination
            out_prefix = Path(OUT_DIR) / sex / tissue / gene / locus
            out_prefix.parent.mkdir(parents=True, exist_ok=True)

            # Write GWAS z-score file
            gwas_z = out_prefix.parent / f"{locus}_gwas.z"
            grp[["rsid", "gwas_zscore"]].drop_duplicates("rsid").to_csv(
                gwas_z, sep="\t", index=False, header=False
            )

            # Write eQTL z-score file
            eqtl_z = out_prefix.parent / f"{locus}_eqtl.z"
            grp[["rsid", "eqtl_zscore"]].drop_duplicates("rsid").to_csv(
                eqtl_z, sep="\t", index=False, header=False
            )

            # Copy/create LD file
            ld_out = out_prefix.parent / f"{locus}.ld"
            if os.path.exists(ld_path):
                import shutil
                shutil.copy(ld_path, ld_out)
            else:
                # 1x1 identity if no LD file
                np.savetxt(ld_out, np.array([[1.0]]), fmt="%.6f", delimiter="\t")

            summary_rows.append({
                "sex": sex, "tissue": tissue, "gene": gene, "locus": locus
            })

summary_df = pd.DataFrame(summary_rows).drop_duplicates()
summary_df.to_csv(SUMMARY_OUT, index=False)
print(f"Prepared {len(summary_df)} locus-tissue-gene combinations")
print(f"Summary saved to: {SUMMARY_OUT}")
print(f"\nBreakdown:")
print(summary_df.groupby(["sex","tissue"]).size().to_string())
