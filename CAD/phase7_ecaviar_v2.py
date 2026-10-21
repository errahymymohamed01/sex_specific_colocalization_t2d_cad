#!/usr/bin/env python3
"""
eCAVIAR Python implementation v2 — CAD version
Handles single-SNP loci (1x1 LD matrix) which is the case for CAD.
Uses Wakefield Approximate Bayes Factors (W=0.15).
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
parser.add_argument("--prior",   type=float, default=0.01)
parser.add_argument("--W",       type=float, default=0.15)
args = parser.parse_args()

Path(args.outdir).mkdir(parents=True, exist_ok=True)


def wakefield_abf(z, W=0.15):
    """Wakefield Approximate Bayes Factor for a single z-score."""
    V = 1.0
    log_abf = 0.5 * np.log(V / (V + W)) + 0.5 * (W / (V + W)) * (z ** 2)
    return np.exp(log_abf)


def wakefield_pip(z_scores, prior=0.01, W=0.15):
    """Compute PIP for multiple SNPs."""
    z = np.array(z_scores, dtype=float)
    V = 1.0
    log_abf = 0.5 * np.log(V / (V + W)) + 0.5 * (W / (V + W)) * (z ** 2 / V)
    log_post = np.log(prior) + log_abf
    log_post -= log_post.max()
    post = np.exp(log_post)
    return post / post.sum()


def run_locus(gwas_z_file, eqtl_z_file, out_prefix, prior, W):
    try:
        gdf = pd.read_csv(gwas_z_file, sep="\t", header=None, names=["rsid", "z"])
        edf = pd.read_csv(eqtl_z_file, sep="\t", header=None, names=["rsid", "z"])
    except Exception as e:
        return None, f"Load error: {e}"

    # Common SNPs
    common = list(dict.fromkeys(s for s in gdf["rsid"] if s in set(edf["rsid"])))

    if len(common) == 0:
        return None, "No common SNPs"

    gz = gdf.set_index("rsid").loc[common, "z"].values.astype(float)
    ez = edf.set_index("rsid").loc[common, "z"].values.astype(float)
    gz = np.nan_to_num(gz, nan=0.0)
    ez = np.nan_to_num(ez, nan=0.0)

    if len(common) == 1:
        # Single-SNP locus: CLPP = PIP_gwas * PIP_eqtl directly from ABFs
        abf_g = wakefield_abf(gz[0], W)
        abf_e = wakefield_abf(ez[0], W)
        pip_g = np.array([(prior * abf_g) / (1 - prior + prior * abf_g)])
        pip_e = np.array([(prior * abf_e) / (1 - prior + prior * abf_e)])
    else:
        pip_g = wakefield_pip(gz, prior, W)
        pip_e = wakefield_pip(ez, prior, W)

    clpp = float(np.sum(pip_g * pip_e))

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


# ── Main ──────────────────────────────────────────────────
summary = pd.read_csv(args.summary)
print(f"Running eCAVIAR (Wakefield ABF) on {len(summary)} combinations")
print(f"Prior={args.prior}  W={args.W}\n")

results = []
failed  = []
done    = 0

for _, row in summary.iterrows():
    sex    = row["sex"]
    locus  = row["locus"]
    tissue = row["tissue"]
    gene   = row["gene"]

    locus_dir  = Path(args.inputs) / sex / tissue / gene
    out_prefix = Path(args.outdir) / sex / tissue / gene / locus

    gwas_z = locus_dir / f"{locus}_gwas.z"
    eqtl_z = locus_dir / f"{locus}_eqtl.z"

    if not (gwas_z.exists() and eqtl_z.exists()):
        failed.append({"sex":sex,"locus":locus,"tissue":tissue,"gene":gene,"error":"missing files"})
        continue

    clpp, err = run_locus(str(gwas_z), str(eqtl_z), str(out_prefix), args.prior, args.W)
    done += 1

    if done % 20 == 0 or done == len(summary):
        clpp_str = f"{clpp:.4f}" if clpp is not None else "NA"
        print(f"  [{done}/{len(summary)}] {sex}/{tissue}/{locus}  CLPP={clpp_str}")

    if err:
        failed.append({"sex":sex,"locus":locus,"tissue":tissue,"gene":gene,"error":err})
    else:
        results.append({
            "sex": sex, "locus": locus, "tissue": tissue, "gene": gene,
            "clpp": clpp,
            "colocalized": clpp >= 0.01
        })

# ── Results ───────────────────────────────────────────────
if not results:
    print("No results — all combinations failed.")
    if failed:
        print("Failures:")
        for f in failed[:10]:
            print(f"  {f}")
else:
    results_df = pd.DataFrame(results).sort_values("clpp", ascending=False)
    results_df.to_csv(Path(args.outdir) / "clpp_all_results.csv", index=False, float_format="%.6f")

    print(f"\n{'='*55}")
    print(f"Completed: {len(results)}/{len(summary)}")
    print(f"Failed:    {len(failed)}")
    print(f"Colocalized (CLPP >= 0.01): {results_df['colocalized'].sum()}")
    print(f"Colocalized (CLPP >= 0.10): {(results_df['clpp'] >= 0.1).sum()}")

    print(f"\nTop 20 results:")
    print(results_df.head(20)[["sex","locus","tissue","gene","clpp"]].to_string(index=False))

    print(f"\nBy sex:")
    print(results_df.groupby("sex")["colocalized"].sum().to_string())

    print(f"\nBy tissue:")
    print(results_df.groupby("tissue")["colocalized"].sum().to_string())

    print(f"\nSaved to: {args.outdir}/clpp_all_results.csv")
