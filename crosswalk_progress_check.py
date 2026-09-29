
import pandas as pd
from pathlib import Path

CROSSWALK_PATH = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\accession_crosswalk.csv")


def is_filled(series: pd.Series) -> pd.Series:
    """Robust check: True where a real (non-missing, non-blank) value exists."""
    return series.notna() & (series.astype(str).str.strip() != "") & \
           (series.astype(str).str.strip().str.lower() != "nan")


def main():
    df = pd.read_csv(CROSSWALK_PATH)

    filled_mask = is_filled(df["ipk_gbis_id"])
    confirmed = df[filled_mask]
    pending = df[~filled_mask]

    print("=" * 80)
    print("CROSSWALK COMPLETION STATUS")
    print("=" * 80)
    print(f"Total accessions: {len(df)}")
    print(f"Confirmed (ipk_gbis_id filled): {len(confirmed)} ({100*len(confirmed)/len(df):.1f}%)")
    print(f"Still pending: {len(pending)} ({100*len(pending)/len(df):.1f}%)")

    print(f"\nBy identifier type:")
    for id_type, group in df.groupby("parsed_id_type"):
        n_filled = is_filled(group["ipk_gbis_id"]).sum()
        print(f"  {id_type}: {n_filled}/{len(group)} confirmed")

    if len(confirmed) > 0:
        dupes = confirmed[confirmed.duplicated(subset=["ipk_gbis_id"], keep=False)]
        if not dupes.empty:
            print(f"\n*** WARNING: {len(dupes)} rows share a duplicate ipk_gbis_id ***")
            print(dupes[["assemblytics_accession", "ipk_gbis_id"]].to_string(index=False))
        else:
            print(f"\nNo duplicate GBIS ID assignments detected.")

    if len(pending) > 0:
        print(f"\nStill pending ({len(pending)}):")
        print(pending[["assemblytics_accession", "parsed_id_type", "gbis_search_url"]].head(20).to_string(index=False))
        if len(pending) > 20:
            print(f"  ... and {len(pending) - 20} more")


if __name__ == "__main__":
    main()