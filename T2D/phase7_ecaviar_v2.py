#!/usr/bin/env python3
"""
eCAVIAR Python implementation v2
Uses Wakefield Approximate Bayes Factors (W=0.15, standard for GWAS)
with proper prior for fine-mapping context (not genome-wide).
"""

import argparse, warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--inputs",  required=True)
parser.add_argument("--outdir",  default="ecaviar_output")
parser.add_argument("--summary", required=True)
parser.add_argument("--prior",   type=float, default=0.01,
                    help="Prior per SNP (default 0.01 for fine-mapping context)")
parser.add_argument("--W",       type=float, default=0.15,
                    help="Wakefield prior variance (default 0.15)")
args = parser.parse_args()

Path(args.outdir).mkdir(parents=True, exist_ok=True)


def wakefield_pip(z_scores, prior=0.01, W=0.15):
    """
    Compute PIP using Wakefield Approximate Bayes Factors.
    
    log ABF_i = 0.5 * log(V/(V+W)) + 0.5 * (W/(V+W)) * z_i²
    where V=1 (standardized z-scores), W=prior effect variance.
    
    PIP_i = prior * ABF_i / (1 - prior + prior * ABF_i)
    
    Then normalize so sum(PIP) = 1 (fine-mapping normalization).
    """
    z = np.array(z_scores, dtype=float)
    V = 1.0

    log_abf = 0.5 * np.log(V / (V + W)) + 0.5 * (W / (V + W)) * (z ** 2 / V)

    # Posterior: proportional to prior * ABF
    log_post = np.log(prior) + log_abf

    # Normalize (log-sum-exp)
    log_post -= log_post.max()
    post = np.exp(log_post)
    pip = post / post.sum()

    return pip


def run_locus(gwas_z_file, eqtl_z_file, ld_file, out_prefix, prior, W):
    try:
        gdf = pd.read_csv(gwas_z_file, sep="\t", header=None, names=["rsid","z"])
        edf = pd.read_csv(eqtl_z_file, sep="\t", header=None, names=["rsid","z"])
    except Exception as e:
        return None, f"Load error: {e}"

    # Common SNPs
    common = list(dict.fromkeys(s for s in gdf["rsid"] if s in set(edf["rsid"])))
    if len(common) < 2:
        return None, f"Only {len(common)} common SNPs"

    gz = gdf.set_index("rsid").loc[common, "z"].values.astype(float)
    ez = edf.set_index("rsid").loc[common, "z"].values.astype(float)

    # Replace NaN
    gz = np.nan_to_num(gz, nan=0.0)
    ez = np.nan_to_num(ez, nan=0.0)

    pip_g = wakefield_pip(gz, prior, W)
    pip_e = wakefield_pip(ez, prior, W)
    clpp  = float(np.sum(pip_g * pip_e))

    # Save per-SNP results
    result_df = pd.DataFrame({
        "rsid":        common,
        "gwas_z":      gz,
        "eqtl_z":      ez,
        "pip_gwas":    pip_g,
        "pip_eqtl":    pip_e,
        "pip_product": pip_g * pip_e,
    }).sort_values("pip_product", ascending=False)

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(f"{out_prefix}_snp.tsv", sep="\t", index=False, float_format="%.6f")

    return clpp, None


# ── Main ──────────────────────────────────────────────────────────────────────
summary = pd.read_csv(args.summary)
print(f"Running eCAVIAR (Wakefield ABF) on {len(summary)} combinations")
print(f"Prior={args.prior}  W={args.W}\n")

results = []
failed  = []
done    = 0

for _, row in summary.iterrows():
    sex, locus, tissue, gene = row["sex"], row["locus"], row["tissue"], row["gene"]

    prefix     = Path(args.inputs) / sex / tissue / gene / locus
    out_prefix = Path(args.outdir) / sex / tissue / gene / locus

    gwas_z = f"{prefix}_gwas.z"
    eqtl_z = f"{prefix}_eqtl.z"
    ld_f   = f"{prefix}.ld"

    if not (Path(gwas_z).exists() and Path(eqtl_z).exists()):
        failed.append({"sex":sex,"locus":locus,"tissue":tissue,"gene":gene,"error":"missing"})
        continue

    clpp, err = run_locus(gwas_z, eqtl_z, ld_f, str(out_prefix), args.prior, args.W)
    done += 1

    if done % 100 == 0 or done == len(summary):
        clpp_str = f"{clpp:.4f}" if clpp is not None else "NA"
        print(f"  [{done}/{len(summary)}] {sex}/{tissue}/{locus}  CLPP={clpp_str}")

    if err:
        failed.append({"sex":sex,"locus":locus,"tissue":tissue,"gene":gene,"error":err})
    else:
        results.append({"sex":sex,"locus":locus,"tissue":tissue,"gene":gene,"clpp":clpp,
                        "colocalized": clpp >= 0.01})

# ── Results ───────────────────────────────────────────────────────────────────
results_df = pd.DataFrame(results).sort_values("clpp", ascending=False)
results_df.to_csv(Path(args.outdir) / "clpp_all_results.csv", index=False, float_format="%.6f")

print(f"\n{'='*55}")
print(f"Completed: {len(results)}/{len(summary)}")
print(f"Colocalized (CLPP >= 0.01): {results_df['colocalized'].sum()}")
print(f"Colocalized (CLPP >= 0.10): {(results_df['clpp'] >= 0.1).sum()}")

print(f"\nTop 20 results:")
print(results_df.head(20)[["sex","locus","tissue","gene","clpp"]].to_string(index=False))

print(f"\nBy sex:")
print(results_df.groupby("sex")["colocalized"].sum().to_string())

print(f"\nBy tissue:")
print(results_df.groupby("tissue")["colocalized"].sum().to_string())

print(f"\nSaved to: {args.outdir}/clpp_all_results.csv")
