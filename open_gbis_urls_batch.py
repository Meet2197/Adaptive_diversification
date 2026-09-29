
import pandas as pd
import webbrowser
import time
from pathlib import Path

CROSSWALK_PATH = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\accession_crosswalk.csv")
BATCH_SIZE = 10  # opens 10 tabs at a time -- manageable to work through


def main():
    df = pd.read_csv(CROSSWALK_PATH)

    # Only open ones not yet confirmed (ipk_gbis_id still blank)
    pending = df[df["ipk_gbis_id"].isna() | (df["ipk_gbis_id"].astype(str).str.strip() == "")]
    print(f"Total accessions: {len(df)}")
    print(f"Still pending GBIS confirmation: {len(pending)}")

    if pending.empty:
        print("All accessions already confirmed!")
        return

    batch = pending.head(BATCH_SIZE)
    print(f"\nOpening {len(batch)} URLs in browser (batch of {BATCH_SIZE}):\n")

    for _, row in batch.iterrows():
        print(f"  {row['assemblytics_accession']:30s} -> {row['gbis_search_url']}")
        webbrowser.open_new_tab(row["gbis_search_url"])
        time.sleep(1.0)  # stagger tab opening

    print(f"\n{len(batch)} tabs opened. After confirming each in GBIS, fill in "
          f"'ipk_gbis_id' in {CROSSWALK_PATH}, save, then re-run this script "
          f"for the next batch.")


if __name__ == "__main__":
    main()