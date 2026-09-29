
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

BASE_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")
QC_REPORT = BASE_DIR / "01_QC" / "dataset_QC_report.csv"
OUTPUT_PATH = BASE_DIR.parent / "accession_crosswalk.csv"  # matches extended_analysis.py's expected location

GBIS_SEARCH_BASE = "https://gbis.ipk-gatersleben.de/gbis2/r/WebFrontend?query="


def parse_accession_name(name: str) -> dict:
    """
    Classifies an Assemblytics accession name (e.g. 'MorexV3_HOR_3474')
    into its underlying genebank identifier type and number/name.
    """
    stripped = re.sub(r"^MorexV3_", "", name, flags=re.IGNORECASE)

    hor_match = re.match(r"^HOR[_\s]?(\d+)$", stripped, flags=re.IGNORECASE)
    if hor_match:
        hor_number = hor_match.group(1)
        return {
            "id_type": "HOR",
            "id_value": hor_number,
            "gbis_query": f"HOR {hor_number}",
        }

    wbdc_match = re.match(r"^WBDC[_\s]?(\d+)$", stripped, flags=re.IGNORECASE)
    if wbdc_match:
        wbdc_number = wbdc_match.group(1)
        return {
            "id_type": "WBDC",
            "id_value": wbdc_number,
            "gbis_query": f"WBDC {wbdc_number}",
        }

    # Fallback: treat remaining underscores as part of a compound cultivar
    # name (e.g. 'Chiba_ZDM02064', 'Chikurin_Ibaraki') -- flagged for
    # manual review rather than guessed at.
    return {
        "id_type": "named_cultivar_or_unrecognized_pattern",
        "id_value": stripped,
        "gbis_query": stripped.replace("_", " "),
    }


def build_crosswalk_template() -> pd.DataFrame:
    if not QC_REPORT.exists():
        raise FileNotFoundError(
            f"{QC_REPORT} not found. Run main.py (Step 1) first -- this "
            f"tool reads the accession list directly from its QC output."
        )

    qc = pd.read_csv(QC_REPORT)
    if "accession" not in qc.columns:
        raise ValueError(f"'accession' column not found in {QC_REPORT}")

    accessions = sorted(qc["accession"].unique().tolist())
    print(f"Found {len(accessions)} accessions in {QC_REPORT}")

    rows = []
    for acc in accessions:
        parsed = parse_accession_name(acc)
        gbis_url = GBIS_SEARCH_BASE + parsed["gbis_query"].replace(" ", "+")
        rows.append({
            "assemblytics_accession": acc,
            "parsed_id_type": parsed["id_type"],
            "parsed_id_value": parsed["id_value"],
            "gbis_search_url": gbis_url,
            "ipk_gbis_id": "",              # fill in manually after checking GBIS
            "microbiome_genotype_id": "",    # fill in manually after checking PRJEB50020/21 metadata
            "notes": "",
        })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = build_crosswalk_template()

    df.to_csv(OUTPUT_PATH, index=False)

    print(f"\nCrosswalk template written to: {OUTPUT_PATH}")
    print(f"Total accessions: {len(df)}")
    print("\nBreakdown by identifier type:")
    print(df["parsed_id_type"].value_counts().to_string())

    unrecognized = df[df["parsed_id_type"] == "named_cultivar_or_unrecognized_pattern"]
    if len(unrecognized) > 0:
        print(f"\n{len(unrecognized)} accessions did not match HOR_#### or WBDC_### "
              f"patterns (likely named cultivars, e.g. 'Bowman', 'Barke', or "
              f"compound names). These need manual GBIS lookup by cultivar "
              f"name instead of accession number -- their gbis_search_url "
              f"has been generated using the name directly.")
        print(unrecognized[["assemblytics_accession", "parsed_id_value"]].to_string(index=False))

    print(
        "\nNEXT STEP: Open each 'gbis_search_url' in a browser, confirm the "
        "matching IPK GBIS record, and manually fill in 'ipk_gbis_id' in "
        f"{OUTPUT_PATH}. Do NOT auto-fill these -- GBIS confirmation "
        "requires visual verification (taxonomy, donor, collection date) "
        "to avoid mismatches between similarly-numbered accessions."
    )