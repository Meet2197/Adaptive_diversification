
import subprocess
from pathlib import Path
import requests

BASE_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")
BAM_DIR = BASE_DIR / "external_data" / "bams"
READS_DIR = BASE_DIR / "external_data" / "reads"
BAM_DIR.mkdir(parents=True, exist_ok=True)
READS_DIR.mkdir(parents=True, exist_ok=True)

ENA_PORTAL_FILEREPORT_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# ------------------------------------------------------------------
# FILL THIS IN once you've confirmed real run accessions from ERZ6810918
# or its linked study. Format: {your_accession_name: ena_run_accession}
# Example (PLACEHOLDER ONLY -- these are NOT verified real accessions):
#   RUN_ACCESSIONS = {
#       "HOR_3474": "ERR1234567",
#       "WBDC_237": "ERR1234568",
#   }
# ------------------------------------------------------------------
RUN_ACCESSIONS: dict[str, str] = {
    # <-- populate this after confirming real ENA run accessions
}


def get_file_links(run_accession: str) -> dict:
    """Query ENA Portal API for FASTQ/BAM/CRAM FTP links for a given run."""
    params = {
        "accession": run_accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,submitted_ftp,submitted_format",
        "format": "tsv",
    }
    response = requests.get(ENA_PORTAL_FILEREPORT_API, params=params, timeout=30)
    response.raise_for_status()

    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        raise ValueError(f"No file report returned for {run_accession}")

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    return dict(zip(header, values))


def download_ftp_file(ftp_url: str, out_path: Path):
    """Download a single file from an ENA FTP link (given as host/path, no scheme)."""
    full_url = f"https://{ftp_url}" if not ftp_url.startswith("http") else ftp_url
    print(f"    Downloading: {full_url}")
    response = requests.get(full_url, stream=True, timeout=120)
    response.raise_for_status()

    with open(out_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"    Saved: {out_path} ({out_path.stat().st_size:,} bytes)")


def sort_and_index_bam(unsorted_bam: Path) -> Path:
    """Sort and index a BAM file using samtools (must be installed / on PATH)."""
    sorted_path = unsorted_bam.with_name(unsorted_bam.stem + ".sorted.bam")
    subprocess.run(["samtools", "sort", "-o", str(sorted_path), str(unsorted_bam)], check=True)
    subprocess.run(["samtools", "index", str(sorted_path)], check=True)
    return sorted_path


def main():
    if not RUN_ACCESSIONS:
        print(
            "RUN_ACCESSIONS is empty. Populate it with confirmed ENA run "
            "accessions before running this script. See the docstring for "
            "instructions on how to find them (via ERZ6810918's linked "
            "study/sample records)."
        )
        return

    print("=" * 80)
    print("FETCHING READ DATA FROM ENA")
    print("=" * 80)

    for local_accession, run_accession in RUN_ACCESSIONS.items():
        print(f"\n{local_accession} -> {run_accession}")
        try:
            file_info = get_file_links(run_accession)
        except (requests.exceptions.RequestException, ValueError) as e:
            print(f"  ERROR fetching file report: {e}")
            continue

        submitted_ftp = file_info.get("submitted_ftp", "")
        submitted_format = file_info.get("submitted_format", "").upper()
        fastq_ftp = file_info.get("fastq_ftp", "")

        if "BAM" in submitted_format or "CRAM" in submitted_format:
            for ftp_link in submitted_ftp.split(";"):
                if not ftp_link:
                    continue
                out_name = f"{local_accession}_{Path(ftp_link).name}"
                out_path = BAM_DIR / out_name
                try:
                    download_ftp_file(ftp_link, out_path)
                except requests.exceptions.RequestException as e:
                    print(f"  ERROR downloading {ftp_link}: {e}")
        elif fastq_ftp:
            print(
                f"  NOTE: only FASTQ available for {run_accession} "
                f"(no pre-aligned BAM/CRAM). Downloading raw reads to "
                f"{READS_DIR} -- these will require alignment to MorexV3 "
                f"before use in Step 10."
            )
            for ftp_link in fastq_ftp.split(";"):
                if not ftp_link:
                    continue
                out_name = f"{local_accession}_{Path(ftp_link).name}"
                out_path = READS_DIR / out_name
                try:
                    download_ftp_file(ftp_link, out_path)
                except requests.exceptions.RequestException as e:
                    print(f"  ERROR downloading {ftp_link}: {e}")
        else:
            print(f"  WARNING: no fastq_ftp or submitted_ftp found for {run_accession}")

    print(
        "\nIf files were downloaded to 'reads/' (FASTQ only), you still "
        "need to align them to MorexV3 (e.g., using bwa-mem2 or minimap2) "
        "and sort/index with samtools before Step 10 read-validation can "
        "use them. That alignment step is not automated here -- confirm "
        "read availability first."
    )


if __name__ == "__main__":
    main()