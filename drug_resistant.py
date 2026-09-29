import pandas as pd
import re
from pathlib import Path

TABLE6_PATH = Path(r"C:\Users\bhatt\Desktop\Omics\diversification_barley\table6.xlsx")
QC_REPORT = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis\01_QC\dataset_QC_report.csv")
CROSSWALK_PATH = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\accession_crosswalk.csv")


def extract_id(name):
    name = re.sub(r"^MorexV3_", "", name, flags=re.IGNORECASE)
    m = re.match(r"^(HOR|WBDC)_?(\d+)$", name, flags=re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}_{m.group(2)}"
    return None


def main():
    full_table = pd.read_excel(TABLE6_PATH, sheet_name="single_plant_per_accession")
    print(f"Full identity table: {len(full_table):,} rows")
    print(f"Columns: {full_table.columns.tolist()}\n")

    qc = pd.read_csv(QC_REPORT)
    your_accessions = qc["accession"].unique().tolist()

    your_df = pd.DataFrame({
        "assemblytics_accession": your_accessions,
        "extracted_id": [extract_id(a) for a in your_accessions],
    })

    merged = your_df.merge(
        full_table, left_on="extracted_id", right_on="accession", how="left", indicator=True
    )

    n_matched = (merged["_merge"] == "both").sum()
    print("=" * 80)
    print("FULL 22,626-ROW MATCH RESULTS")
    print("=" * 80)
    print(f"Matched: {n_matched} / {len(your_df)}")
    print(f"Not matched: {len(your_df) - n_matched} / {len(your_df)}\n")

    matched = merged[merged["_merge"] == "both"]
    print(matched[["assemblytics_accession", "accession", "panel", "name",
                    "country_of_origin", "status_genebank_information"]].to_string(index=False))

    unmatched = merged[merged["_merge"] != "both"]
    print(f"\n--- Still not matched ({len(unmatched)}) ---")
    print(unmatched[["assemblytics_accession", "extracted_id"]].to_string(index=False))

    # Update crosswalk with ALL confirmed matches (union of Table 1 + this table)
    crosswalk = pd.read_csv(CROSSWALK_PATH)
    crosswalk["ipk_gbis_id"] = crosswalk["ipk_gbis_id"].astype("object")

    update_map = dict(zip(matched["assemblytics_accession"], matched["accession"]))
    n_updated = 0
    for idx, row in crosswalk.iterrows():
        current = row.get("ipk_gbis_id")
        is_blank = pd.isna(current) or str(current).strip() == ""
        if is_blank and row["assemblytics_accession"] in update_map:
            crosswalk.at[idx, "ipk_gbis_id"] = update_map[row["assemblytics_accession"]]
            n_updated += 1

    crosswalk.to_csv(CROSSWALK_PATH, index=False)
    print(f"\n{n_updated} NEW rows updated (in addition to prior 16 from Table 1)")
    print(f"Total confirmed now: {(crosswalk['ipk_gbis_id'].notna() & (crosswalk['ipk_gbis_id'].astype(str).str.strip() != '')).sum()} / 74")


if __name__ == "__main__":
    main()