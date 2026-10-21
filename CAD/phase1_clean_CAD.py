import os
import csv
import gzip

TRAIT      = "CAD"
INPUT_DIR  = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD"
OUTPUT_DIR = r"C:\Users\HP\Desktop\master BIAM\s2\BioInformatics I\Projet\CAD\Cleaned"

KEEP_COLS = ["rsid", "CHR", "BP", "reference_allele", "NEA",
             "eaf", "beta", "se", "z_score", "p_value", "n"]

def clean_file(input_path, output_path):
    kept = removed = no_rsid = no_allele = 0

    with open(input_path, 'r', encoding='utf-8') as f_in, \
         gzip.open(output_path, 'wt', encoding='utf-8') as f_out:

        reader = csv.DictReader(f_in)
        writer = csv.DictWriter(f_out, fieldnames=KEEP_COLS)
        writer.writeheader()

        for row in reader:
            rsid = row.get("rsid", "").strip()
            if not rsid or rsid == ".":
                no_rsid += 1; removed += 1; continue

            ea = row.get("reference_allele", "").strip()
            if not ea:
                no_allele += 1; removed += 1; continue

            try:
                beta = float(row["beta"])
                se   = float(row["se"])
                if se == 0:
                    removed += 1; continue
            except:
                removed += 1; continue

            z_score = round(beta / se, 4)

            writer.writerow({
                "rsid":             rsid,
                "CHR":              row.get("CHR", "").strip(),
                "BP":               row.get("BP", "").strip(),
                "reference_allele": ea,
                "NEA":              ".",
                "eaf":              row.get("eaf", "").strip(),
                "beta":             beta,
                "se":               se,
                "z_score":          z_score,
                "p_value":          row.get("p_value", "").strip(),
                "n":                row.get("n", "").strip()
            })
            kept += 1

    return kept, removed, no_rsid, no_allele

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_files = [f for f in os.listdir(INPUT_DIR)
                 if f.startswith(TRAIT + "_") and f.endswith(".tsv")]

    if not all_files:
        print(f"No files found for trait: {TRAIT}")
        return

    print(f"\n{'='*60}")
    print(f"  PHASE 1 — Cleaning trait: {TRAIT}")
    print(f"  Files found: {len(all_files)}")
    print(f"  NOTE: NEA set to '.' — will be added from dbSNP later")
    print(f"{'='*60}\n")

    for filename in sorted(all_files):
        input_path  = os.path.join(INPUT_DIR, filename)
        output_path = os.path.join(OUTPUT_DIR,
                      filename.replace(".tsv", "_clean.tsv.gz"))

        print(f"  Processing: {filename}")
        kept, removed, no_rsid, no_allele = clean_file(input_path, output_path)

        print(f"     Kept:    {kept:,}")
        print(f"     Removed: {removed:,}  (no rsid: {no_rsid} | no allele: {no_allele})")
        print(f"     Saved:   {output_path}\n")

    print(f"{'='*60}")
    print(f"  DONE — {TRAIT} Phase 1 complete!")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
