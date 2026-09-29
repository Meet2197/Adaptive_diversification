
import time
from pathlib import Path
import requests

BASE_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")
MICROBIOME_DIR = BASE_DIR / "14_Microbiome_Integration"
XML_DIR = MICROBIOME_DIR / "xml"
MICROBIOME_DIR.mkdir(parents=True, exist_ok=True)
XML_DIR.mkdir(parents=True, exist_ok=True)

# ENA public REST APIs (documented, stable endpoints):
ENA_BROWSER_XML_API = "https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
ENA_PORTAL_FILEREPORT_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# PRJEB50061 (parent) + confirmed component sub-projects from the ENA page
PROJECT_ACCESSIONS = {
    "PRJEB50061": "Parent project - Identifying plant genes shaping microbiota composition in the barley rhizosphere",
    "PRJEB50019": "16S rRNA, bulk soil + rhizosphere, NAM HEB-25 Family 15 (QTL mapping)",
    "PRJEB50020": "16S rRNA, elite + sibling lines (QTL validation) -- PRIORITY: most likely genotype overlap",
    "PRJEB50021": "ITS rRNA (fungal), elite + sibling lines (QTL validation)",
    "PRJEB50400": "Root RNA-seq, controlled conditions (transcriptome, not microbiome composition)",
}


def fetch_project_xml(accession: str) -> Path:
    """Fetch raw XML metadata for a project/study accession via ENA Browser API."""
    url = ENA_BROWSER_XML_API.format(accession=accession)
    print(f"  Fetching XML: {url}")
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    out_path = XML_DIR / f"{accession}.xml"
    out_path.write_bytes(response.content)
    print(f"    Saved: {out_path} ({len(response.content):,} bytes)")
    return out_path


def fetch_sample_report_tsv(accession: str) -> Path:
    """
    Fetch a sample/run-level TSV report via ENA Portal API. Requests common
    fields useful for building the accession crosswalk (genotype/sample
    alias, run accession, FASTQ FTP links).
    """
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,sample_accession,sample_alias,"
                  "sample_title,scientific_name,library_strategy,"
                  "fastq_ftp,submitted_ftp,study_accession",
        "format": "tsv",
    }
    print(f"  Fetching sample report for: {accession}")
    response = requests.get(ENA_PORTAL_FILEREPORT_API, params=params, timeout=30)
    response.raise_for_status()

    out_path = MICROBIOME_DIR / f"{accession}_sample_report.tsv"
    out_path.write_text(response.text, encoding="utf-8")

    n_lines = response.text.count("\n")
    print(f"    Saved: {out_path} ({n_lines} lines, incl. header)")
    return out_path


def main():
    print("=" * 80)
    print("FETCHING PRJEB50061 MICROBIOME PROJECT METADATA (ENA)")
    print("=" * 80)

    results = []
    for accession, description in PROJECT_ACCESSIONS.items():
        print(f"\n{accession}: {description}")
        try:
            xml_path = fetch_project_xml(accession)
            tsv_path = fetch_sample_report_tsv(accession)
            results.append({
                "accession": accession, "status": "success",
                "xml_path": str(xml_path), "tsv_path": str(tsv_path),
            })
        except requests.exceptions.RequestException as e:
            print(f"    ERROR: {e}")
            results.append({
                "accession": accession, "status": "failed", "error": str(e),
            })
        time.sleep(1)  # be polite to the API

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for r in results:
        print(f"  {r['accession']}: {r['status']}")

    print(
        "\nNEXT STEPS:\n"
        "  1. Open each *_sample_report.tsv in "
        f"{MICROBIOME_DIR} and inspect the 'sample_alias'/'sample_title' "
        "columns for genotype names.\n"
        "  2. Cross-check those genotype names against your 74 "
        "Assemblytics accession names (paste both lists to me for "
        "automated comparison).\n"
        "  3. The actual OTU/ASV abundance table is NOT retrievable via "
        "this script -- download it manually from the Supplementary "
        "Information of Escudero-Martinez et al. 2022 "
        "(https://doi.org/10.1038/s41467-022-30115-y) and place it at:\n"
        f"     {MICROBIOME_DIR / 'PRJEB50061_otu_table.csv'}\n"
    )


if __name__ == "__main__":
    main()