
from __future__ import annotations
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
import logging

log = logging.getLogger("ena_xml_tools")


def flatten_xml(path: Path) -> pd.DataFrame:
    
    tree = ET.parse(path)
    root = tree.getroot()

    rows = []

    def _walk(elem, parent_tag="ROOT"):
        text = (elem.text or "").strip()
        rows.append({
            "tag": elem.tag,
            "parent_tag": parent_tag,
            "attributes": dict(elem.attrib),
            "text": text if text else None,
        })
        for child in elem:
            _walk(child, parent_tag=elem.tag)

    _walk(root)
    return pd.DataFrame(rows)


def find_file_references(flat_df: pd.DataFrame) -> pd.DataFrame:
    candidates = []
    for _, row in flat_df.iterrows():
        attrs = row["attributes"]
        tag_lower = row["tag"].lower()

        is_file_like_tag = any(
            key in tag_lower for key in ["file", "url", "link"]
        )
        has_file_like_attr = any(
            any(k in str(key).lower() for k in ["filename", "url", "ftp", "checksum", "filetype"])
            for key in attrs.keys()
        )

        if is_file_like_tag or has_file_like_attr:
            candidates.append({
                "tag": row["tag"],
                "parent_tag": row["parent_tag"],
                "attributes": attrs,
                "text": row["text"],
            })

    return pd.DataFrame(candidates)


def find_accession_references(flat_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    accession_like_keys = {"accession", "alias", "refname"}

    for _, row in flat_df.iterrows():
        attrs = row["attributes"]
        for key, value in attrs.items():
            if str(key).lower() in accession_like_keys:
                rows.append({
                    "tag": row["tag"],
                    "attribute_key": key,
                    "value": value,
                })

    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def summarize_ena_xml(path: Path, label: str, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    flat = flatten_xml(path)
    flat_out = output_dir / f"{label}_flattened_structure.csv"
    flat.to_csv(flat_out, index=False)

    file_refs = find_file_references(flat)
    file_refs_out = output_dir / f"{label}_file_references.csv"
    file_refs.to_csv(file_refs_out, index=False)

    acc_refs = find_accession_references(flat)
    acc_refs_out = output_dir / f"{label}_accession_references.csv"
    acc_refs.to_csv(acc_refs_out, index=False)

    summary = {
        "source_xml": str(path),
        "n_total_elements": len(flat),
        "n_file_like_references": len(file_refs),
        "n_accession_references": len(acc_refs),
        "unique_tags_present": sorted(flat["tag"].unique().tolist()),
        "outputs": {
            "flattened_structure": str(flat_out),
            "file_references": str(file_refs_out),
            "accession_references": str(acc_refs_out),
        },
    }

    log.info(f"[{label}] Parsed {path.name}: {summary['n_total_elements']} elements, "
              f"{summary['n_file_like_references']} file-like references, "
              f"{summary['n_accession_references']} accession references.")

    if summary["n_file_like_references"] == 0:
        log.warning(
            f"[{label}] No file/FTP/URL references found in this XML. "
            f"This record likely describes METADATA ONLY (e.g. a project "
            f"or study summary) -- actual sequence/data files are NOT "
            f"embedded here and would need separate retrieval from ENA "
            f"(e.g. via the project's run/sample accessions)."
        )
    return summary