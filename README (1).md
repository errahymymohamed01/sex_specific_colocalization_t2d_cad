# sex_specific_colocalization_t2d_cad

Sex-stratified colocalization of GWAS and eQTL signals for Type 2 Diabetes (T2D) and Coronary Artery Disease (CAD), using a custom Python re-implementation of eCAVIAR based on Wakefield's Approximate Bayes Factors (Master BIAM, Bioinformatics I, 2025-2026).

# Sex-Specific Genetic Colocalization in T2D and CAD

Most GWAS pool men and women into one cohort, which can hide variants that act mainly in one sex. This project runs the same 8-phase pipeline on **female-only** and **male-only** GWAS summary statistics, and links them to **GTEx v8 eQTLs** in 7 relevant tissues to find which genes and tissues may drive each sex's disease risk.

## Overview

| Trait | Role in the project | LD information |
|-------|---------------------|----------------|
| T2D   | Primary analysis    | Real 200x200 LD matrices (LDlink, 1000G CEU) |
| CAD   | Secondary, preliminary analysis | 1x1 identity matrices (LDlink `ldmatrix` server failure) |

**Authors:** Mohamed Errahymy, Chaimae Rami
**Supervisor:** Prof. Cherng Yang

## Data

- **GWAS summary statistics (GRCh37/hg19), split by sex**
  - T2D: about 7,000,000 SNPs per sex
  - CAD: 7,426,763 SNPs per file (both / female / male)
  - Fields used: `rsid`, `CHR`, `BP`, effect allele, NEA, `eaf`, `beta`, `se`, `z_score`, `p_value`, `n`
  - The CAD raw file had no NEA column, which requires an extra recovery step (Phase 2)
- **GTEx v8 eQTLs** (positional IDs `CHR_BP_REF_ALT_b37`, with `gtex_beta`, `gtex_se`, `gtex_pval`)
  - Adipose Subcutaneous, Adipose Visceral (Omentum), Artery Aorta, Artery Coronary, Artery Tibial, Liver, Pancreas
  - The same 7 tissues are used for both traits
- **Reference panels:** dbSNP b151 GRCh37 (NEA annotation), 1000 Genomes CEU via the NIH LDlink API (LD)

## Pipeline

| Phase | Step | What it does | T2D script | CAD script |
|-------|------|--------------|------------|------------|
| 1 | Cleaning | Standardizes chromosome names (RefSeq to integers), drops missing values and `se = 0`, computes `z = beta / se`, splits by sex | `phase1_clean.py` | `phase1_clean_CAD.py` |
| 2 | NEA annotation (CAD only) | Recovers the missing non-effect allele from dbSNP b151 (all ALT alleles + complement-strand fallback), run on Google Colab | not needed | `phase2_Untitled4.ipynb` |
| 3 | Locus selection | Distance-based clumping: lead SNP at p < 5e-8, +/- 500 kb exclusion window, processed per chromosome / in chunks | `phase3_clump_T2D.py` | `phase3_clump_CAD.py` |
| 4 | SNP matching | Matches GWAS rsIDs to GTEx positional IDs, aligns alleles, flips z-score sign when needed, drops palindromic (A/T, C/G) and mismatched SNPs | `phase4_rematch_loci.py` | `phase4_snp_match_CAD.py` |
| 5 | LD matrices | Pairwise R2 from LDlink (LDproxy + LDmatrix, batches of 50 SNPs, retry with backoff, resumable). Alternative: local 1000G VCF + PLINK | `phase5ecaviar_ld_pipeline.py` | `phase5_ld_CAD.py` |
| 6 | eCAVIAR inputs | Builds synchronized `.gwas.z`, `.eqtl.z` and `.ld` files per locus / tissue / gene, with strict SNP ordering | `phase6_prepare_ecaviar_inputs.py` | `phase6_prepare_ecaviar_CAD.py` |
| 7 | Colocalization | Wakefield ABF, then PIP, then CLPP | `phase7_ecaviar_v2.py` | `phase7_ecaviar_v2.py` |
| 8 | Sex comparison | Classifies each locus-gene pair as female-only, male-only or shared | `phase8_analyze_clpp.py` (in `CAD/`) | `phase8_analyze_clpp.py` |

### Statistical core (Phase 7)

```
log ABF_i = 0.5 * log(V / (V + W)) + 0.5 * (W / (V + W)) * z_i^2 / V
PIP_i     = (pi0 * ABF_i) / sum_j (pi0 * ABF_j)
CLPP      = sum_i PIP_i(GWAS) * PIP_i(eQTL)
```

Parameters: `V = 1.0`, `W = 0.15` (Wakefield 2009), fine-mapping prior `pi0 = 0.01`. A standard genome-wide prior (1e-4) pushed every CLPP to about 0, so the prior was set for a fine-mapping context (inside a window already known to contain a significant signal).

## Technical challenges handled

- **Memory:** chunked reading (100,000 rows) of 650 MB GWAS files, keeping only SNPs inside locus windows (over 99% reduction)
- **LDlink API:** token must be a URL query parameter (in the POST body it fails silently with HTTP 200), batches over 50 SNPs trigger HTTP 500
- **ID and format mismatches:** rsIDs vs positional IDs, X/Y/MT chromosome strings in GTEx, no precomputed eQTL z-scores
- **Allele direction:** sign correction of z-scores when GWAS and GTEx alleles are swapped
- **LD matrix order:** SNP order in the z-score files must match the LD matrix exactly

## Key Results

**T2D (fully validated, real LD)**

- 75 female and 102 male independent loci, 1,226 locus-tissue-gene combinations tested
- 1,039 unique colocalized pairs (CLPP >= 0.01): 53.6% female-only, 45.5% male-only, 0.9% shared
- Strongest signals in adipose tissue; top hit rs3132520 / LINC01697 in Adipose Visceral (CLPP = 0.410)

**CAD (preliminary, 1x1 LD workaround)**

- 12 female and 61 male independent loci, 117 combinations tested
- 68 unique colocalized pairs: 10.3% female-only, 77.9% male-only, 11.8% shared
- Notable leads: rs12740374 / SORT1 in Liver, rs9349379 / PHACTR1 in Artery Tibial
- Without real LD, CAD CLPP scores are upper-bound approximations, not confirmed colocalizations

## Limitations

- CAD fine-mapping has no real LD information (LDlink `ldmatrix` returned HTTP 500 for every CAD locus)
- CAD final sample is small (68 pairs), so percentages can shift with a few loci
- LD computed on European (CEU) 1000G only
- 200-SNP-per-locus cap may truncate large regulatory regions
- Single causal variant assumption (Wakefield ABF); multi-signal loci would need FINEMAP or similar
- The CLPP values depend on the prior (pi0 = 0.01)
- GTEx is bulk-tissue and not sex-stratified

## Repository Structure

```
T2D/
  phase1_clean.py
  phase3_clump_T2D.py
  phase4_rematch_loci.py
  phase5ecaviar_ld_pipeline.py
  phase6_prepare_ecaviar_inputs.py
  phase7_ecaviar_v2.py

CAD/
  phase1_clean_CAD.py
  phase2_Untitled4.ipynb            (NEA recovery, Google Colab)
  phase3_clump_CAD.py
  phase4_snp_match_CAD.py
  phase5_ld_CAD.py
  phase6_prepare_ecaviar_CAD.py
  phase7_ecaviar_v2.py
  phase8_analyze_clpp.py

report/                              full LaTeX report (PDF)
```

Notes:
- T2D has no Phase 2 because its raw file already contains the NEA column.
- Phase 8 (sex comparison) is applied to the merged CLPP results of both traits.
- Both folders use the same `phase7_ecaviar_v2.py` logic (Wakefield ABF, `pi0 = 0.01`).

## How to Run

1. Place the cleaned sex-stratified GWAS files and the GTEx v8 eQTL files in a `data/` folder, and update the paths at the top of each script.
2. Get a free LDlink API token (NIH) and set it as an environment variable.
3. Install the dependencies:

```
pip install numpy pandas requests matplotlib seaborn scipy
```

4. Run the scripts of one trait in phase order:

```
# T2D
cd T2D
python phase1_clean.py
python phase3_clump_T2D.py
python phase4_rematch_loci.py
python phase5ecaviar_ld_pipeline.py
python phase6_prepare_ecaviar_inputs.py
python phase7_ecaviar_v2.py

# CAD
cd CAD
python phase1_clean_CAD.py
# phase2_Untitled4.ipynb  -> run on Google Colab (needs the 15 GB dbSNP VCF)
python phase3_clump_CAD.py
python phase4_snp_match_CAD.py
python phase5_ld_CAD.py
python phase6_prepare_ecaviar_CAD.py
python phase7_ecaviar_v2.py
python phase8_analyze_clpp.py
```

Phase 5 is resumable: already generated LD matrices are skipped if the run is interrupted.

## Next Steps

- Recompute CAD LD locally (PLINK `--r2` on 1000G EUR) to replace the 1x1 matrices
- Validate top CAD hits with independent eQTL datasets
- Extend to other cardiometabolic traits (lipids, blood pressure, obesity)
- Cross-sex validation and single-cell / trans-eQTL resolution

## References

- Hormozdiari et al. (2016). Colocalization of GWAS and eQTL signals detects target genes. *AJHG*
- Wakefield (2009). Bayes factors for genome-wide association studies. *Genetic Epidemiology*
- GTEx Consortium (2020). *Science*
- Machiela & Chanock (2015). LDlink. *Bioinformatics*
- Mahajan et al. (2022). Fine-mapping global variation in type 2 diabetes susceptibility. *Nature Genetics*

## Author

Mohamed Errahymy, Master BIAM
[GitHub](https://github.com/errahymymohamed01)
