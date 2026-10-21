import gzip
import os
import csv

# ============================================================
# CONFIGURATION
# ============================================================

TRAIT      = "T2D"
INPUT_DIR  = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\T2D"
OUTPUT_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\T2D\Cleaned"

KEEP_COLS = ["rsid", "CHR", "BP", "reference_allele", "NEA", "eaf", "beta", "se", "z_score", "p_value", "n"]

# ============================================================
# CHROM converter: NC_000001.10 → 1
# ============================================================
def convert_chrom(chrom_str):
    mapping = {
        "NC_000001.10": "1",  "NC_000002.11": "2",  "NC_000003.11": "3",
        "NC_000004.11": "4",  "NC_000005.09": "5",  "NC_000006.11": "6",
        "NC_000007.13": "7",  "NC_000008.10": "8",  "NC_000009.11": "9",
        "NC_000010.10": "10", "NC_000011.09": "11", "NC_000012.11": "12",
        "NC_000013.10": "13", "NC_000014.08": "14", "NC_000015.09": "15",
        "NC_000016.09": "16", "NC_000017.10": "17", "NC_000018.09": "18",
        "NC_000019.09": "19", "NC_000020.10": "20", "NC_000021.08": "21",
        "NC_000022.10": "22", "NC_000023.10": "X",  "NC_000024.09": "Y"
    }
    return mapping.get(chrom_str, chrom_str)

# ============================================================
# CLEAN ONE FILE
# ============================================================
def clean_file(input_path, output_path):
    kept      = 0
    removed   = 0
    no_rsid   = 0
    no_allele = 0

    with gzip.open(input_path, 'rt', encoding='utf-8') as f_in, \
         gzip.open(output_path, 'wt', encoding='utf-8') as f_out:

        reader = csv.DictReader(f_in)
        writer = csv.DictWriter(f_out, fieldnames=KEEP_COLS)
        writer.writeheader()

        for row in reader:

            rsid = row.get("rsid", "").strip()
            if not rsid or rsid == ".":
                no_rsid += 1
                removed += 1
                continue

            ea  = row.get("reference_allele", "").strip()
            nea = row.get("NEA", "").strip()
            if not ea or not nea:
                no_allele += 1
                removed += 1
                continue

            try:
                beta = float(row["beta"])
                se   = float(row["se"])
                if se == 0:
                    removed += 1
                    continue
            except:
                removed += 1
                continue

            z_score     = round(beta / se, 4)
            chrom_fixed = convert_chrom(row.get("CHROM", "").strip())

            writer.writerow({
                "rsid":             rsid,
                "CHR":              row.get("CHR", chrom_fixed).strip(),
                "BP":               row.get("BP", "").strip(),
                "reference_allele": ea,
                "NEA":              nea,
                "eaf":              row.get("eaf", "").strip(),
                "beta":             beta,
                "se":               se,
                "z_score":          z_score,
                "p_value":          row.get("p_value", "").strip(),
                "n":                row.get("n", "").strip()
            })
            kept += 1

    return kept, removed, no_rsid, no_allele

# ============================================================
# MAIN
# ============================================================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_files = [f for f in os.listdir(INPUT_DIR) if f.startswith(TRAIT + "_") and f.endswith(".tsv.gz")]

    if not all_files:
        print(f"❌ No files found for trait: {TRAIT}")
        return

    print(f"\n{'='*60}")
    print(f"  PHASE 1 — Cleaning trait: {TRAIT}")
    print(f"  Files found: {len(all_files)}")
    print(f"{'='*60}\n")

    for filename in sorted(all_files):
        input_path  = os.path.join(INPUT_DIR, filename)
        output_path = os.path.join(OUTPUT_DIR, filename.replace(".tsv.gz", "_clean.tsv.gz"))

        print(f"  ⏳ Processing: {filename}")
        kept, removed, no_rsid, no_allele = clean_file(input_path, output_path)

        print(f"     ✅ Kept:    {kept:,}")
        print(f"     ❌ Removed: {removed:,}  (no rsid: {no_rsid} | no allele: {no_allele})")
        print(f"     💾 Saved to: {output_path}\n")

    print(f"{'='*60}")
    print(f"  ✅ DONE — {TRAIT} cleaning complete!")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()