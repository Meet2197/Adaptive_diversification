from __future__ import annotations
from ena_xml_tools import summarize_ena_xml
import argparse
import json
import logging
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from concurrent.futures import ProcessPoolExecutor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger("extended_analysis")


# ============================================================
# 0. CONFIGURATION
# ============================================================

BASE_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")

DIRS = {
    "final": BASE_DIR / "06_Final_Datasets",
    "catalogue": BASE_DIR / "07_SV_Catalogue",
    "genotype": BASE_DIR / "08_Genotype_Matrix",
    "phenotype": BASE_DIR / "09_Phenotype_Integration",
    "gwas": BASE_DIR / "10_PanGWAS_Results",
    "annotation": BASE_DIR / "11_Functional_Annotation",
    "graph": BASE_DIR / "12_Graph_Pangenome",
    "validation": BASE_DIR / "13_Read_Validation",
    "microbiome": BASE_DIR / "14_Microbiome_Integration",
    "manifest": BASE_DIR / "00_Manifest",
    "external": BASE_DIR / "external_data",
}
for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)

TOOLS = {
    "survivor": shutil.which("SURVIVOR") or "SURVIVOR",
    "gemma": shutil.which("gemma") or "gemma",
    "minigraph": shutil.which("minigraph") or "minigraph",
    "cactus": shutil.which("cactus-pangenome") or "cactus-pangenome",
    "vg": shutil.which("vg") or "vg",
    "samtools": shutil.which("samtools") or "samtools",
}

REQUIRED_INPUT_FILES = {
    "cleaned_sv": DIRS["final"] / "ALL_Assemblytics_SVs_CLEANED.csv",
    "accession_crosswalk": BASE_DIR / "accession_crosswalk.csv",
    "phenotype_table": DIRS["external"] / "drought_phenotypes.csv",
    "gff3_annotation": DIRS["external"] / "genes.gff3",
    "microbiome_otu_table": DIRS["external"] / "PRJEB50061_otu_table.csv",
    "microbiome_metadata": DIRS["external"] / "PRJEB50061_sample_metadata.csv",
    "bam_dir": DIRS["external"] / "bams",
    "bam_metadata_xml": DIRS["external"] / "5. BAM_files_ERZ6810918.xml",
    "microbiome_metadata_xml": DIRS["external"] / "6. Microbiome_PRJEB50061.xml",
}
KNOWN_PLACED_CHROMOSOMES = {f"chr{i}H" for i in range(1, 8)}


def require_file(key: str) -> Path:
    path = REQUIRED_INPUT_FILES[key]
    if not path.exists():
        raise FileNotFoundError(
            f"Required input '{key}' not found at: {path}\n"
            f"This step cannot produce a real result without it."
        )
    return path


# ============================================================
# 1. MANIFEST LOGGER
# ============================================================

@dataclass
class ManifestEntry:
    step: str
    description: str
    inputs: list = field(default_factory=list)
    outputs: list = field(default_factory=list)
    parameters: dict = field(default_factory=dict)
    status: str = "not_run"
    software_versions: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    notes: str = ""


class Manifest:
    def __init__(self, path: Path = DIRS["manifest"] / "analysis_manifest.jsonl"):
        self.path = path

    def log(self, entry: ManifestEntry):
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry.__dict__) + "\n")
        log.info(f"[MANIFEST] {entry.step}: {entry.status}")


MANIFEST = Manifest()


# ============================================================
# STEP 3a — SCALABLE SV CLUSTERING (connected components)
# ============================================================

def _cluster_one_group(
    idx: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    sizes: np.ndarray,
    is_insertion: bool,
    max_distance_bp: int,
    min_reciprocal_overlap: float,
    max_size_ratio_diff: float,
) -> tuple:
    order = np.argsort(starts, kind="mergesort")
    idx, starts, ends, sizes = idx[order], starts[order], ends[order], sizes[order]
    n = len(idx)
    if n < 2:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    right_bounds = np.searchsorted(starts, starts + max_distance_bp, side="right")
    edge_i_list, edge_j_list = [], []

    for i in range(n - 1):
        j_end = right_bounds[i]
        if j_end <= i + 1:
            continue
        js = np.arange(i + 1, j_end)

        size_i = max(sizes[i], 1)
        size_j = np.maximum(sizes[js], 1)
        size_ratio_diff = np.abs(sizes[js] - sizes[i]) / np.maximum(size_i, size_j)

        if is_insertion:
            match = size_ratio_diff <= max_size_ratio_diff
        else:
            inter = np.minimum(ends[js], ends[i]) - np.maximum(starts[js], starts[i])
            inter = np.clip(inter, 0, None)
            len_i = max(ends[i] - starts[i], 1)
            len_j = np.maximum(ends[js] - starts[js], 1)
            reciprocal_overlap = np.minimum(inter / len_i, inter / len_j)
            match = (reciprocal_overlap >= min_reciprocal_overlap) & (
                size_ratio_diff <= max_size_ratio_diff
            )

        matched_js = js[match]
        if len(matched_js):
            edge_i_list.extend([idx[i]] * len(matched_js))
            edge_j_list.extend(idx[matched_js].tolist())

    return (
        np.array(edge_i_list, dtype=np.int64),
        np.array(edge_j_list, dtype=np.int64),
    )


def cluster_svs_scalable(
    df: pd.DataFrame,
    max_distance_bp: int = 500,
    min_reciprocal_overlap: float = 0.7,
    max_size_ratio_diff: float = 0.3,
    chrom_col: str = "reference",
    start_col: str = "ref_start",
    end_col: str = "ref_stop",
    type_col: str = "SV_TYPE_NORMALIZED",
    size_col: str = "SV_SIZE",
    n_workers: int = 1,
) -> pd.DataFrame:
    df = df.reset_index(drop=True)
    groups = list(df.groupby([chrom_col, type_col], sort=False))
    log.info(f"Clustering {len(df):,} records across {len(groups)} (chrom,type) groups")

    def _run(group_key_and_df):
        (chrom, sv_type), gdf = group_key_and_df
        is_insertion = str(sv_type).lower() == "insertion"
        return _cluster_one_group(
            gdf.index.to_numpy(),
            gdf[start_col].to_numpy(),
            gdf[end_col].to_numpy(),
            gdf[size_col].to_numpy(),
            is_insertion,
            max_distance_bp,
            min_reciprocal_overlap,
            max_size_ratio_diff,
        )

    all_edges_i, all_edges_j = [], []

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for ei, ej in ex.map(_run, groups):
                all_edges_i.append(ei)
                all_edges_j.append(ej)
    else:
        for i, g in enumerate(groups):
            ei, ej = _run(g)
            all_edges_i.append(ei)
            all_edges_j.append(ej)
            if (i + 1) % 20 == 0:
                log.info(f"  processed {i+1}/{len(groups)} groups")

    edge_i = np.concatenate(all_edges_i) if all_edges_i else np.array([], dtype=np.int64)
    edge_j = np.concatenate(all_edges_j) if all_edges_j else np.array([], dtype=np.int64)

    n_total = len(df)
    graph = coo_matrix(
        (np.ones(len(edge_i), dtype=np.int8), (edge_i, edge_j)),
        shape=(n_total, n_total),
    )
    n_components, labels = connected_components(csgraph=graph, directed=False)
    log.info(f"Raw connected components (PRE chain-split fix): {n_components:,}")

    df = df.copy()
    df["CLUSTER_ID"] = labels
    return df


# ============================================================
# STEP 3b — CHAIN-SPLIT FIX (REQUIRED — closes single-linkage
# chaining loophole confirmed on real chr1H data: 25.5% of
# records were in clusters spanning >1500bp with genuine
# multi-record chaining)
# ============================================================

def split_chained_clusters(
    df: pd.DataFrame,
    max_distance_bp: int = 500,
    span_flag_multiplier: float = 3.0,
    max_plausible_sv_size: int = 1_000_000,
) -> pd.DataFrame:
    """
    Re-splits clusters affected by single-linkage chaining, so that within
    any resulting sub-cluster, max(ref_start) - min(ref_start) <=
    max_distance_bp (equivalent guarantee to complete-linkage clustering
    with a hard distance cap, for 1-D breakpoint coordinates).

    PERFORMANCE NOTE: this implementation avoids repeated O(n_total) scans
    of the full dataframe per cluster (the original bottleneck). Cluster
    membership is looked up via a precomputed index dict (O(1) per
    cluster), and new cluster IDs are tracked incrementally instead of
    recomputing df["CLUSTER_ID"].max() on every iteration.
    """
    df = df.reset_index(drop=True)

    cluster_group = df.groupby("CLUSTER_ID")
    cluster_stats = cluster_group.agg(
        n_records=("SV_RECORD_ID", "count"),
        n_accessions=("accession", "nunique"),
        min_start=("ref_start", "min"),
        max_end=("ref_stop", "max"),
    )
    cluster_stats["span_bp"] = cluster_stats["max_end"] - cluster_stats["min_start"]
    cluster_stats["records_per_accession"] = (
        cluster_stats["n_records"] / cluster_stats["n_accessions"]
    )

    oversized_singleton_ids = set(cluster_stats[
        (cluster_stats["n_records"] == 1) &
        (cluster_stats["span_bp"] > max_plausible_sv_size)
    ].index)
    df["FLAGGED_OVERSIZED_SINGLE_SV"] = df["CLUSTER_ID"].isin(oversized_singleton_ids)
    if len(oversized_singleton_ids) > 0:
        log.warning(
            f"{len(oversized_singleton_ids)} single-record SV(s) exceed "
            f"{max_plausible_sv_size}bp reference span. Flagged in "
            f"'FLAGGED_OVERSIZED_SINGLE_SV' -- NOT auto-corrected."
        )

    span_threshold = max_distance_bp * span_flag_multiplier
    flagged_ids = cluster_stats[
        (cluster_stats["n_records"] > 1) &
        (
            (cluster_stats["span_bp"] > span_threshold) |
            (cluster_stats["records_per_accession"] > 1.3)
        )
    ].index

    log.info(
        f"Flagged {len(flagged_ids)} of {len(cluster_stats)} MULTI-RECORD "
        f"clusters for re-splitting (span > {span_threshold}bp or "
        f"same-accession multiplicity > 1.3)"
    )

    n_total = len(df)
    chain_split_flag = np.zeros(n_total, dtype=bool)

    if len(flagged_ids) == 0:
        df["CHAIN_SPLIT_APPLIED"] = chain_split_flag
        return df

    # --- KEY FIX 1: build cluster -> row-index lookup ONCE (O(n_total)),
    #     instead of re-scanning the full column on every iteration ---
    group_indices = cluster_group.indices  # dict: cluster_id -> np.ndarray of row positions

    cluster_id_arr = df["CLUSTER_ID"].to_numpy().copy()
    start_arr = df["ref_start"].to_numpy()

    # --- KEY FIX 2: track next_id incrementally, never rescan full array ---
    next_id = int(cluster_id_arr.max()) + 1

    t_start = time.time()
    n_done = 0

    for cid in flagged_ids:
        idx = group_indices[cid]
        if len(idx) < 2:
            continue

        vals = start_arr[idx]
        order = np.argsort(vals, kind="mergesort")
        sorted_idx = idx[order]
        sorted_vals = vals[order]

        # --- KEY FIX 3: simple O(k) greedy 1-D diameter-constrained split
        #     instead of scipy pdist/linkage/fcluster per cluster. For 1-D
        #     data this gives the same guarantee (no two members in a
        #     resulting sub-cluster exceed max_distance_bp apart), with
        #     far less per-call overhead. ---
        sub_labels = np.empty(len(sorted_vals), dtype=np.int64)
        current_label = 0
        group_start_val = sorted_vals[0]
        sub_labels[0] = current_label
        for k in range(1, len(sorted_vals)):
            if sorted_vals[k] - group_start_val > max_distance_bp:
                current_label += 1
                group_start_val = sorted_vals[k]
            sub_labels[k] = current_label

        new_ids = next_id + sub_labels
        cluster_id_arr[sorted_idx] = new_ids
        chain_split_flag[sorted_idx] = True
        next_id = int(new_ids.max()) + 1

        n_done += 1
        if n_done % 20000 == 0:
            elapsed = time.time() - t_start
            rate = n_done / elapsed
            eta_sec = (len(flagged_ids) - n_done) / rate if rate > 0 else float("nan")
            log.info(
                f"  chain-split progress: {n_done}/{len(flagged_ids)} clusters "
                f"processed ({rate:.0f}/sec, ETA {eta_sec/60:.1f} min)"
            )

    df["CLUSTER_ID"] = cluster_id_arr
    df["CHAIN_SPLIT_APPLIED"] = chain_split_flag
    return df


def compare_clustering_thresholds_scalable(
    df: pd.DataFrame,
    distance_grid=(100, 250, 500, 1000, 2000),
    overlap_grid=(0.5, 0.7, 0.9),
    n_workers: int = 1,
) -> pd.DataFrame:
    rows = []
    n_accessions = df["accession"].nunique()
    for dist in distance_grid:
        for ro in overlap_grid:
            t0 = time.time()
            clustered = cluster_svs_scalable(
                df, max_distance_bp=dist, min_reciprocal_overlap=ro, n_workers=n_workers
            )
            clustered = split_chained_clusters(clustered, max_distance_bp=dist)
            support = clustered.groupby("CLUSTER_ID")["accession"].nunique()
            row = {
                "max_distance_bp": dist,
                "min_reciprocal_overlap": ro,
                "n_clusters": clustered["CLUSTER_ID"].nunique(),
                "n_private": int((support == 1).sum()),
                "n_shared": int((support > 1).sum()),
                "n_core_all_accessions": int((support == n_accessions).sum()),
                "n_chain_split_records": int(clustered["CHAIN_SPLIT_APPLIED"].sum()),
                "runtime_sec": round(time.time() - t0, 1),
            }
            log.info(row)
            rows.append(row)
    return pd.DataFrame(rows)


# ============================================================
# STEP 4 — NON-REDUNDANT SV CATALOGUE
# ============================================================

def build_sv_catalogue(clustered_df: pd.DataFrame) -> pd.DataFrame:
    """
    One row per CLUSTER_ID (non-redundant SV locus). Frequency classes are
    defined relative to the ACTUAL evaluated accession set, not the barley
    species as a whole:
        private                       -> supported by exactly 1 accession
        shared                        -> supported by >1 but not all
        core_all_evaluated_accessions -> supported by every accession
    """
    n_total_accessions = clustered_df["accession"].nunique()

    has_chain_cols = (
        "CHAIN_SPLIT_APPLIED" in clustered_df.columns and
        "FLAGGED_OVERSIZED_SINGLE_SV" in clustered_df.columns
    )

    agg_kwargs = dict(
        reference=("reference", "first"),
        start=("ref_start", "min"),
        end=("ref_stop", "max"),
        SV_type=("SV_TYPE_NORMALIZED", "first"),
        median_size=("SV_SIZE", "median"),
        n_supporting_accessions=("accession", "nunique"),
        supporting_accessions=("accession", lambda s: ";".join(sorted(set(s)))),
        source_records=("SV_RECORD_ID", lambda s: ";".join(s)),
    )
    if has_chain_cols:
        agg_kwargs["was_chain_split"] = ("CHAIN_SPLIT_APPLIED", "any")
        agg_kwargs["contains_oversized_flag"] = ("FLAGGED_OVERSIZED_SINGLE_SV", "any")

    agg = (
        clustered_df.groupby("CLUSTER_ID")
        .agg(**agg_kwargs)
        .reset_index()
        .sort_values("CLUSTER_ID")
        .reset_index(drop=True)
    )

    agg["SV_ID"] = [f"SVLOC_{i+1:07d}" for i in range(len(agg))]

    def classify(n):
        if n == 1:
            return "private"
        if n == n_total_accessions:
            return "core_all_evaluated_accessions"
        return "shared"

    agg["FREQUENCY_CLASS"] = agg["n_supporting_accessions"].apply(classify)
    agg["n_total_accessions_evaluated"] = n_total_accessions
    coordinate_semantics_map = {
        "Insertion": "tight_breakpoint_interval",
        "Deletion": "tight_breakpoint_interval",
        "Tandem_expansion": "tight_breakpoint_interval",
        "Tandem_contraction": "tight_breakpoint_interval",
        "Repeat_expansion": "ambiguous_alignment_window_size_is_copy_number_estimate",
        "Repeat_contraction": "ambiguous_alignment_window_size_is_copy_number_estimate",
    }
    agg["COORDINATE_SEMANTICS"] = agg["SV_type"].map(coordinate_semantics_map).fillna("unknown")
    agg["IS_PLACED_CHROMOSOME"] = agg["reference"].isin(KNOWN_PLACED_CHROMOSOMES)

    output_cols = [
        "SV_ID", "CLUSTER_ID", "reference", "start", "end", "SV_type", "median_size",
        "n_supporting_accessions", "n_total_accessions_evaluated",
        "FREQUENCY_CLASS", "COORDINATE_SEMANTICS", "IS_PLACED_CHROMOSOME",
    ]
    if has_chain_cols:
        output_cols += ["was_chain_split", "contains_oversized_flag"]
    output_cols += ["supporting_accessions", "source_records"]

    return agg[output_cols]


# ============================================================
# STEP 5 — ACCESSION x SV GENOTYPE MATRIX
# ============================================================

def build_genotype_matrix(catalogue: pd.DataFrame, clustered_df: pd.DataFrame):
    accessions = sorted(clustered_df["accession"].unique())

    presence = (
        clustered_df.groupby(["CLUSTER_ID", "accession"]).size().unstack(fill_value=0)
    )
    presence = presence.clip(upper=1)
    presence = presence.reindex(columns=accessions, fill_value=0)

    cluster_to_svid = catalogue.set_index("CLUSTER_ID")["SV_ID"]
    presence.index = presence.index.map(cluster_to_svid)
    presence.index.name = "SV_ID"

    unfiltered = presence.astype(float)

    call_rate = unfiltered.notna().mean(axis=1)
    alt_freq = unfiltered.mean(axis=1, skipna=True)
    mvf = np.minimum(alt_freq, 1 - alt_freq)

    quality_filtered = unfiltered.loc[(call_rate >= 0.8) & (mvf > 0)]

    qc_table = pd.DataFrame({
        "SV_ID": unfiltered.index,
        "call_rate": call_rate.values,
        "alt_frequency": alt_freq.values,
        "minor_variant_frequency": mvf.values,
    })

    return unfiltered, quality_filtered, qc_table


# ============================================================
# STEP 6 — PHENOTYPE INTEGRATION
# ============================================================

def load_phenotype_table() -> pd.DataFrame:
    path = require_file("phenotype_table")
    df = pd.read_csv(path)
    required_cols = {"accession_id", "trait", "value", "unit", "condition",
                      "replicate", "environment"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Phenotype table missing required columns: {missing}.")
    return df


def harmonize_accessions(genotype_ids: list) -> pd.DataFrame:
    crosswalk_path = require_file("accession_crosswalk")
    crosswalk = pd.read_csv(crosswalk_path)
    required = {"assemblytics_accession", "ipk_gbis_id", "microbiome_genotype_id"}
    missing = required - set(crosswalk.columns)
    if missing:
        raise ValueError(f"Crosswalk file missing columns: {missing}")

    overlap = crosswalk[crosswalk["assemblytics_accession"].isin(genotype_ids)]
    report = {
        "n_genotype_accessions": len(genotype_ids),
        "n_crosswalk_entries": len(crosswalk),
        "n_overlap_with_genotypes": len(overlap),
    }
    log.info(f"Accession overlap report: {report}")
    return overlap


def summarize_phenotype_replicates(pheno_df: pd.DataFrame, trait: str) -> pd.DataFrame:
    import statsmodels.formula.api as smf

    sub = pheno_df[pheno_df["trait"] == trait].copy()
    if sub.empty:
        raise ValueError(f"No phenotype records found for trait '{trait}'.")

    formula = "value ~ C(environment) + C(condition)" if sub["environment"].nunique() > 1 \
        else "value ~ C(condition)"

    model = smf.mixedlm(formula, sub, groups=sub["accession_id"])
    fit = model.fit(reml=True)
    random_effects = fit.random_effects
    summary = pd.DataFrame(
        [{"accession_id": k, "trait": trait, "adjusted_value": float(v.iloc[0])}
         for k, v in random_effects.items()]
    )
    return summary


# ============================================================
# STEP 7 — SV-PanGWAS
# ============================================================

def write_bimbam_genotype(genotype_matrix: pd.DataFrame, out_path: Path):
    with open(out_path, "w") as f:
        for sv_id, row in genotype_matrix.iterrows():
            values = ",".join("NA" if pd.isna(v) else str(int(v)) for v in row)
            f.write(f"{sv_id},A,T,{values}\n")


def write_bimbam_phenotype(pheno_series: pd.Series, out_path: Path):
    pheno_series.to_csv(out_path, header=False, index=False, na_rep="NA")


def compute_kinship_matrix(genotype_matrix: pd.DataFrame) -> pd.DataFrame:
    fill_value = genotype_matrix.mean(axis=1, skipna=True).mean()
    X = genotype_matrix.fillna(fill_value).to_numpy().T
    X = X - X.mean(axis=0)
    K = (X @ X.T) / X.shape[1]
    accessions = genotype_matrix.columns.tolist()
    return pd.DataFrame(K, index=accessions, columns=accessions)


def run_gemma_lmm(genotype_bimbam, phenotype_file, kinship_file, output_prefix, outdir=DIRS["gwas"]) -> Path:
    if not shutil.which(TOOLS["gemma"]):
        raise FileNotFoundError("GEMMA executable not found.")
    cmd = [
        TOOLS["gemma"], "-g", str(genotype_bimbam), "-p", str(phenotype_file),
        "-k", str(kinship_file), "-lmm", "4", "-outdir", str(outdir), "-o", output_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"GEMMA failed:\n{result.stderr}")
    return outdir / f"{output_prefix}.assoc.txt"


def parse_gemma_output(assoc_path: Path) -> pd.DataFrame:
    from statsmodels.stats.multitest import multipletests
    df = pd.read_csv(assoc_path, sep="\t")
    _, p_adj, _, _ = multipletests(df["p_lrt"], method="fdr_bh")
    df["p_adj_fdr_bh"] = p_adj
    return df


# ============================================================
# STEP 8 — FUNCTIONAL ANNOTATION
# ============================================================
def harmonize_gff3_chromosome_names(gff3_path: Path, output_path: Path,
                                      keep_featuretypes: set = None) -> Path:
    if output_path.exists():
        log.info(f"Harmonized/filtered GFF3 already exists at {output_path}, skipping re-generation.")
        return output_path

    if keep_featuretypes is None:
        keep_featuretypes = {"gene", "exon"}

    placed_chrom_map = {f"{i}H": f"chr{i}H" for i in range(1, 8)}
    n_lines_written = 0
    n_lines_filtered_out = 0
    n_seqid_remapped = 0

    with open(gff3_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(output_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] in placed_chrom_map:
                        parts[1] = placed_chrom_map[parts[1]]
                        line = " ".join(parts) + "\n"
                fout.write(line)
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                continue

            featuretype = fields[2]
            if featuretype not in keep_featuretypes:
                n_lines_filtered_out += 1
                continue

            if fields[0] in placed_chrom_map:
                fields[0] = placed_chrom_map[fields[0]]
                line = "\t".join(fields)
                n_seqid_remapped += 1

            fout.write(line)
            n_lines_written += 1

    log.info(
        f"Harmonized+filtered GFF3 written to {output_path}: "
        f"{n_lines_written:,} lines kept (featuretypes={keep_featuretypes}), "
        f"{n_lines_filtered_out:,} lines filtered out (other featuretypes), "
        f"{n_seqid_remapped:,} seqid values remapped."
    )
    return output_path


def build_gff_database(force: bool = False):
    import gffutils
    gff3 = require_file("gff3_annotation")

    harmonized_gff3 = DIRS["annotation"] / "MorexV3_annotation_gene_exon_only.gff3"
    harmonize_gff3_chromosome_names(gff3, harmonized_gff3, keep_featuretypes={"gene", "exon"})

    db_path = DIRS["annotation"] / "morexv3_annotation.db"
    if db_path.exists() and not force:
        return gffutils.FeatureDB(str(db_path))

    log.info("Building gffutils database (filtered to gene+exon only, "
              "merge_strategy='create_unique' to avoid the slow duplicate-ID "
              "merge path that caused the previous hang)...")

    return gffutils.create_db(
        str(harmonized_gff3), str(db_path), force=True,
        keep_order=False,
        merge_strategy="create_unique",
        sort_attribute_values=False,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        verbose=True,
        pragmas={"synchronous": "OFF", "journal_mode": "MEMORY"},
    )

def _load_chromosome_feature_arrays(db, chrom: str, featuretype: str):
    starts, ends, strands, ids = [], [], [], []
    for feature in db.region(seqid=chrom, featuretype=featuretype):
        starts.append(feature.start)
        ends.append(feature.end)
        strands.append(feature.strand)
        ids.append(feature.id)

    if not starts:
        return None

    order = np.argsort(starts)
    return {
        "start": np.array(starts)[order],
        "end": np.array(ends)[order],
        "strand": np.array(strands)[order],
        "id": np.array(ids, dtype=object)[order],
    }


def _query_overlaps(arrays: dict, query_start: int, query_end: int) -> list:

    if arrays is None or len(arrays["start"]) == 0:
        return []

    upper = np.searchsorted(arrays["start"], query_end, side="right")
    if upper == 0:
        return []

    candidate_ends = arrays["end"][:upper]
    candidate_starts = arrays["start"][:upper]
    mask = candidate_ends >= query_start

    if not mask.any():
        return []

    return arrays["id"][:upper][mask].tolist()


def _query_promoter_proximity(gene_arrays: dict, query_start: int, query_end: int,
                               promoter_window_bp: int) -> list:
    if gene_arrays is None or len(gene_arrays["start"]) == 0:
        return []

    window_lower = query_start - promoter_window_bp
    window_upper = query_end + promoter_window_bp

    lo = np.searchsorted(gene_arrays["start"], window_lower, side="left")
    hi = np.searchsorted(gene_arrays["start"], window_upper, side="right")
    if lo >= hi:
        return []

    hits = []
    for i in range(lo, hi):
        g_start, g_end, g_strand, g_id = (
            gene_arrays["start"][i], gene_arrays["end"][i],
            gene_arrays["strand"][i], gene_arrays["id"][i],
        )
        if g_strand == "+" and 0 < g_start - query_end <= promoter_window_bp:
            hits.append(g_id)
        elif g_strand == "-" and 0 < query_start - g_end <= promoter_window_bp:
            hits.append(g_id)

    return hits

def build_gff_database(force: bool = False):
    import gffutils
    gff3 = require_file("gff3_annotation")

    harmonized_gff3 = DIRS["annotation"] / "MorexV3_annotation_harmonized.gff3"
    harmonize_gff3_chromosome_names(gff3, harmonized_gff3)

    db_path = DIRS["annotation"] / "morexv3_annotation.db"
    if db_path.exists() and not force:
        return gffutils.FeatureDB(str(db_path))
    return gffutils.create_db(
        str(harmonized_gff3), str(db_path), force=True,
        keep_order=True, merge_strategy="merge", sort_attribute_values=True,
)

def annotate_svs(catalogue: pd.DataFrame, promoter_window_bp: int = 2000) -> pd.DataFrame:
    db = build_gff_database()

    placed_chroms = sorted(catalogue.loc[catalogue["IS_PLACED_CHROMOSOME"], "reference"].unique())
    log.info(f"Pre-loading gene/exon feature arrays for {len(placed_chroms)} placed chromosomes...")

    gene_arrays_by_chrom = {}
    exon_arrays_by_chrom = {}
    for chrom in placed_chroms:
        gene_arrays_by_chrom[chrom] = _load_chromosome_feature_arrays(db, chrom, "gene")
        exon_arrays_by_chrom[chrom] = _load_chromosome_feature_arrays(db, chrom, "exon")
        n_genes = len(gene_arrays_by_chrom[chrom]["start"]) if gene_arrays_by_chrom[chrom] else 0
        n_exons = len(exon_arrays_by_chrom[chrom]["start"]) if exon_arrays_by_chrom[chrom] else 0
        log.info(f"  {chrom}: {n_genes:,} genes, {n_exons:,} exons loaded")

    records = []
    n_total = len(catalogue)
    t0 = time.time()

    for i, sv in enumerate(catalogue.itertuples(index=False), 1):
        if not getattr(sv, "IS_PLACED_CHROMOSOME", True):
            records.append({
                "SV_ID": sv.SV_ID, "gene_overlap": "", "exon_overlap": "",
                "promoter_proximity_genes": "", "evidence_type": "unplaced_scaffold_not_annotated",
            })
            continue

        chrom, start, end = sv.reference, int(sv.start), int(sv.end)
        gene_arrays = gene_arrays_by_chrom.get(chrom)
        exon_arrays = exon_arrays_by_chrom.get(chrom)

        overlaps_gene = _query_overlaps(gene_arrays, start, end)
        overlaps_exon = _query_overlaps(exon_arrays, start, end)
        promoter_hit = _query_promoter_proximity(gene_arrays, start, end, promoter_window_bp)

        # Exclude genes already counted as direct overlaps from promoter list
        promoter_hit = [g for g in promoter_hit if g not in set(overlaps_gene)]

        records.append({
            "SV_ID": sv.SV_ID,
            "gene_overlap": ";".join(sorted(set(overlaps_gene))),
            "exon_overlap": ";".join(sorted(set(overlaps_exon))),
            "promoter_proximity_genes": ";".join(sorted(set(promoter_hit))),
            "evidence_type": "positional_overlap_only",
        })

        if i % 200_000 == 0:
            elapsed = time.time() - t0
            rate = i / elapsed
            eta_min = (n_total - i) / rate / 60
            log.info(f"  annotated {i:,}/{n_total:,} SVs ({rate:.0f}/sec, ETA {eta_min:.1f} min)")

    log.info(f"Annotation complete: {n_total:,} SVs processed in {time.time()-t0:.1f}s")
    return pd.DataFrame(records)


# ============================================================
# STEP 9 — GRAPH-PANGENOME PROJECTION
# ============================================================

def build_minigraph_cactus_graph(assembly_fastas, reference_name, output_prefix, outdir=DIRS["graph"]) -> Path:
    if not shutil.which(TOOLS["cactus"]):
        raise FileNotFoundError("cactus-pangenome not found. Requires Linux/WSL2/Docker.")
    seqfile = outdir / f"{output_prefix}.seqfile.txt"
    with open(seqfile, "w") as f:
        for fasta in assembly_fastas:
            f.write(f"{fasta.stem}\t{fasta}\n")
    cmd = [
        TOOLS["cactus"], str(outdir / f"{output_prefix}_jobstore"), str(seqfile),
        str(outdir / output_prefix), "--reference", reference_name, "--outDir", str(outdir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Minigraph-Cactus failed:\n{result.stderr}")
    return outdir / f"{output_prefix}.gfa"


def project_svs_to_graph(graph_gfa: Path, candidate_bed: Path, outdir: Path = DIRS["graph"]) -> Path:
    if not shutil.which(TOOLS["vg"]):
        raise FileNotFoundError("vg toolkit not found on PATH.")
    vg_path = outdir / (graph_gfa.stem + ".vg")
    with open(vg_path, "wb") as out:
        subprocess.run([TOOLS["vg"], "convert", "-g", str(graph_gfa), "-p"], stdout=out, check=True)
    out_positions = outdir / "candidate_sv_graph_positions.tsv"
    result = subprocess.run(
        [TOOLS["vg"], "find", "-x", str(vg_path), "-E", "-B", str(candidate_bed)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"vg find failed:\n{result.stderr}")
    out_positions.write_text(result.stdout)
    return out_positions


# ============================================================
# STEP 10 — READ-SUPPORT VALIDATION
# ============================================================

def validate_sv_with_reads(bam_path, chrom, start, end, sv_type, min_mapq=20, flank_bp=200) -> dict:
    import pysam
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    split_support = discordant_support = spanning_support = 0
    depths = []
    for read in bam.fetch(chrom, max(0, start - flank_bp), end + flank_bp):
        if read.mapping_quality < min_mapq or read.is_unmapped:
            continue
        depths.append(1)
        if read.has_tag("SA"):
            split_support += 1
        if read.is_paired and not read.is_proper_pair:
            discordant_support += 1
        if read.reference_start <= start and read.reference_end is not None and read.reference_end >= end:
            spanning_support += 1
    bam.close()

    return {
        "chrom": chrom, "start": start, "end": end, "sv_type": sv_type,
        "mean_local_depth": float(np.mean(depths)) if depths else 0.0,
        "split_read_support": split_support,
        "discordant_read_support": discordant_support,
        "spanning_read_support": spanning_support,
        "validation_status": "supported" if (split_support > 0 or spanning_support > 0) else "unsupported_or_untested",
        "note": "Assembly-derived SV != independently validated SV.",
    }


def batch_validate_candidates(candidate_df: pd.DataFrame, bam_dir: Path) -> pd.DataFrame:
    rows = []
    for _, sv in candidate_df.iterrows():
        accession_bam = bam_dir / f"{sv.get('accession', '')}.sorted.bam"
        if not accession_bam.exists():
            rows.append({"SV_ID": sv["SV_ID"], "validation_status": "no_reads_available"})
            continue
        result = validate_sv_with_reads(
            accession_bam, sv["reference"], int(sv["start"]), int(sv["end"]), sv["SV_type"]
        )
        result["SV_ID"] = sv["SV_ID"]
        rows.append(result)
    return pd.DataFrame(rows)


# ============================================================
# SECTION 5 — MICROBIOME INTEGRATION
# ============================================================

def load_microbiome_data():
    otu = pd.read_csv(require_file("microbiome_otu_table"), index_col=0)
    meta = pd.read_csv(require_file("microbiome_metadata"))
    return otu, meta


def compute_alpha_diversity(otu_table: pd.DataFrame) -> pd.Series:
    p = otu_table.div(otu_table.sum(axis=1), axis=0).replace(0, np.nan)
    shannon = -(p * np.log(p)).sum(axis=1, skipna=True)
    return shannon.rename("shannon_diversity")


def clr_transform(otu_table: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    comp = otu_table.replace(0, pseudocount)
    log_comp = np.log(comp)
    return log_comp.sub(log_comp.mean(axis=1), axis=0)


def sv_microbiome_association(genotype_matrix, microbiome_metric, crosswalk, min_overlap_n=15) -> pd.DataFrame:
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    overlap_ids = crosswalk.dropna(subset=["microbiome_genotype_id"])["assemblytics_accession"]
    if len(overlap_ids) < min_overlap_n:
        raise ValueError(
            f"Only {len(overlap_ids)} overlapping accessions (minimum: {min_overlap_n}). "
            f"Treat microbiome study as independent supporting evidence instead."
        )

    y = microbiome_metric.reindex(overlap_ids).dropna()
    results = []
    for sv_id, row in genotype_matrix.loc[:, genotype_matrix.columns.intersection(overlap_ids)].iterrows():
        x = row.reindex(y.index)
        valid = x.notna() & y.notna()
        if valid.sum() < min_overlap_n or x[valid].nunique() < 2:
            continue
        X = sm.add_constant(x[valid].astype(float))
        model = sm.OLS(y[valid].astype(float), X).fit()
        results.append({
            "SV_ID": sv_id, "beta": model.params.get(x.name, np.nan),
            "se": model.bse.get(x.name, np.nan), "p_value": model.pvalues.get(x.name, np.nan),
            "n": int(valid.sum()),
        })
    df = pd.DataFrame(results)
    if not df.empty:
        _, df["p_adj_fdr_bh"], _, _ = multipletests(df["p_value"], method="fdr_bh")
    return df


# ============================================================
# CLI DISPATCH
# ============================================================

def main():
    if len(sys.argv) == 1:
        print(
            "\nERROR: no step specified.\n"
            "Example: python extended_analysis.py cluster --max-distance 500 --min-overlap 0.7 --chrom chr1H\n"
            "Steps: cluster, sensitivity-clustering, catalogue, genotype-matrix, "
            "phenotype-integration, gwas, annotate, graph-projection, read-validation, microbiome\n"
        )
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Extended barley pangenomics pipeline")
    parser.add_argument(
        "step",
        choices=[
            "cluster", "sensitivity-clustering", "catalogue", "genotype-matrix",
            "phenotype-integration", "gwas", "annotate", "graph-projection",
            "read-validation", "microbiome",
            "parse-bam-metadata", "parse-microbiome-metadata",
        ],
    )
    parser.add_argument("--max-distance", type=int, default=500)
    parser.add_argument("--min-overlap", type=float, default=0.7)
    parser.add_argument("--chrom", type=str, default=None)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    if args.step == "cluster":
        sv_path = require_file("cleaned_sv")
        df = pd.read_csv(sv_path)
        if args.chrom:
            df = df[df["reference"].str.contains(args.chrom, case=False, na=False)]
            log.info(f"Restricted to chromosome filter '{args.chrom}': {len(df):,} records")

        t0 = time.time()
        clustered = cluster_svs_scalable(
            df, max_distance_bp=args.max_distance,
            min_reciprocal_overlap=args.min_overlap, n_workers=args.workers,
        )
        n_before = int(clustered["CLUSTER_ID"].nunique())

        clustered = split_chained_clusters(clustered, max_distance_bp=args.max_distance)
        n_after = int(clustered["CLUSTER_ID"].nunique())
        n_chain_split_records = int(clustered["CHAIN_SPLIT_APPLIED"].sum())
        n_oversized_flagged = int(clustered["FLAGGED_OVERSIZED_SINGLE_SV"].sum())

        runtime = round(time.time() - t0, 1)

        suffix = f"_{args.chrom}" if args.chrom else ""
        out = DIRS["catalogue"] / f"clustered_SVs{suffix}.csv"
        clustered.to_csv(out, index=False)

        MANIFEST.log(ManifestEntry(
            step="cluster_svs",
            description="Scalable clustering + single-linkage chain-split fix",
            inputs=[str(sv_path)], outputs=[str(out)],
            parameters={
                "max_distance_bp": args.max_distance,
                "min_reciprocal_overlap": args.min_overlap,
                "chrom_filter": args.chrom,
                "n_workers": args.workers,
                "n_input_records": len(df),
                "n_clusters_before_chain_split": n_before,
                "n_clusters_after_chain_split": n_after,
                "n_records_chain_split_applied": n_chain_split_records,
                "n_oversized_singleton_svs_flagged": n_oversized_flagged,
                "runtime_sec": runtime,
            },
            status="confirmatory",
        ))
        log.info(
            f"Done in {runtime}s. Clusters before split: {n_before:,}, "
            f"after split: {n_after:,}. Chain-split records: {n_chain_split_records:,}. "
            f"Oversized singletons flagged: {n_oversized_flagged:,}."
        )

    elif args.step == "sensitivity-clustering":
        sv_path = require_file("cleaned_sv")
        df = pd.read_csv(sv_path)
        if args.chrom:
            df = df[df["reference"].str.contains(args.chrom, case=False, na=False)]
        else:
            log.warning(
                "Running sensitivity-clustering GENOME-WIDE runs the full "
                "cluster+chain-split pipeline 15 times (~5+ hours at this "
                "dataset's scale, based on the ~20min/genome-wide-run "
                "extrapolation from chr1H timing). Strongly recommend using "
                "--chrom chr1H (or another single chromosome) as a "
                "representative sample for threshold selection instead."
            )
        sens = compare_clustering_thresholds_scalable(df, n_workers=args.workers)
        out = DIRS["catalogue"] / "clustering_sensitivity.csv"
        sens.to_csv(out, index=False)
        MANIFEST.log(ManifestEntry(
            step="sensitivity_clustering", description="Clustering-threshold sensitivity",
            inputs=[str(sv_path)], outputs=[str(out)], status="exploratory",
        ))

    elif args.step == "catalogue":
        clustered_path = DIRS["catalogue"] / "clustered_SVs.csv"
        if not clustered_path.exists():
            raise FileNotFoundError(
                f"{clustered_path} not found. Run genome-wide 'cluster' (no --chrom) first."
            )
        clustered = pd.read_csv(clustered_path)
        if "CHAIN_SPLIT_APPLIED" not in clustered.columns:
            raise RuntimeError(
                "clustered_SVs.csv was produced WITHOUT the chain-split fix "
                "(missing 'CHAIN_SPLIT_APPLIED' column). Re-run 'cluster' with "
                "the current version of this script before building the catalogue."
            )
        catalogue = build_sv_catalogue(clustered)
        out = DIRS["catalogue"] / "nonredundant_SV_catalogue.csv"
        catalogue.to_csv(out, index=False)
        MANIFEST.log(ManifestEntry(
            step="build_catalogue", description="Non-redundant SV catalogue",
            inputs=[str(clustered_path)], outputs=[str(out)],
            parameters={"n_loci": len(catalogue)}, status="confirmatory",
        ))
        log.info(f"Catalogue built: {len(catalogue):,} non-redundant SV loci")

    elif args.step == "genotype-matrix":
        clustered = pd.read_csv(DIRS["catalogue"] / "clustered_SVs.csv")
        catalogue = pd.read_csv(DIRS["catalogue"] / "nonredundant_SV_catalogue.csv")
        unfiltered, filtered, qc = build_genotype_matrix(catalogue, clustered)

        log.info(
            "NOTE: call_rate is currently 1.0 for all loci by construction "
            "(Assemblytics reports only detected SVs; no true no-calls exist "
            "yet). The quality_filtered matrix is therefore currently driven "
            "entirely by minor_variant_frequency > 0, not by call_rate. This "
            "will change once Step 10 read-based validation is available."
        )

        unfiltered.to_csv(DIRS["genotype"] / "genotype_matrix_unfiltered.csv")
        filtered.to_csv(DIRS["genotype"] / "genotype_matrix_quality_filtered.csv")
        qc.to_csv(DIRS["genotype"] / "genotype_matrix_QC.csv", index=False)
        MANIFEST.log(ManifestEntry(
            step="genotype_matrix", description="Accession x SV presence/absence matrix",
            outputs=[str(DIRS["genotype"])],
            parameters={"n_loci_unfiltered": unfiltered.shape[0],
                        "n_loci_quality_filtered": filtered.shape[0],
                        "n_accessions": unfiltered.shape[1]},
            status="confirmatory",
            notes="call_rate=1.0 for all loci (assembly-comparison-derived; no true no-calls yet).",
        ))
        log.info(f"Genotype matrix: {unfiltered.shape[0]:,} loci (unfiltered), "
                  f"{filtered.shape[0]:,} (quality-filtered), {unfiltered.shape[1]} accessions")

    elif args.step == "phenotype-integration":
        genotype_ids = pd.read_csv(
            DIRS["genotype"] / "genotype_matrix_quality_filtered.csv", index_col=0
        ).columns.tolist()
        crosswalk = harmonize_accessions(genotype_ids)
        crosswalk.to_csv(DIRS["phenotype"] / "accession_overlap_report.csv", index=False)
        pheno = load_phenotype_table()
        pheno.to_csv(DIRS["phenotype"] / "phenotype_table_loaded.csv", index=False)
        MANIFEST.log(ManifestEntry(
            step="phenotype_integration", description="Phenotype-genotype harmonization",
            outputs=[str(DIRS["phenotype"])], status="exploratory",
        ))

    elif args.step == "gwas":
        geno = pd.read_csv(DIRS["genotype"] / "genotype_matrix_quality_filtered.csv", index_col=0)
        kinship = compute_kinship_matrix(geno)
        kinship.to_csv(DIRS["gwas"] / "kinship_matrix.csv")
        write_bimbam_genotype(geno, DIRS["gwas"] / "genotype.bimbam")
        MANIFEST.log(ManifestEntry(
            step="gwas_preparation", description="GEMMA input preparation",
            outputs=[str(DIRS["gwas"])], status="exploratory",
            notes="Awaiting harmonized phenotype file; n=74 must be disclosed as effective sample size.",
        ))

    elif args.step == "annotate":
        catalogue = pd.read_csv(DIRS["catalogue"] / "nonredundant_SV_catalogue.csv")
        annotated = annotate_svs(catalogue)
        out = DIRS["annotation"] / "SV_functional_annotation.csv"
        annotated.to_csv(out, index=False)
        MANIFEST.log(ManifestEntry(
            step="annotate_svs", description="GFF3-based positional annotation",
            outputs=[str(out)], status="confirmatory",
            notes="Positional overlap only; not experimentally demonstrated functional effect.",
        ))

    elif args.step == "read-validation":
        catalogue = pd.read_csv(DIRS["catalogue"] / "nonredundant_SV_catalogue.csv")
        bam_dir = require_file("bam_dir")
        result = batch_validate_candidates(catalogue, bam_dir)
        out = DIRS["validation"] / "read_support_validation.csv"
        result.to_csv(out, index=False)
        MANIFEST.log(ManifestEntry(
            step="read_validation", description="Read-based SV validation",
            outputs=[str(out)], status="confirmatory",
        ))

    elif args.step == "microbiome":
        otu, meta = load_microbiome_data()
        alpha = compute_alpha_diversity(otu)
        alpha.to_csv(DIRS["microbiome"] / "alpha_diversity.csv")
        MANIFEST.log(ManifestEntry(
            step="microbiome_diversity", description="Alpha diversity (PRJEB50061)",
            outputs=[str(DIRS["microbiome"])], status="exploratory",
        ))
    elif args.step == "parse-bam-metadata":
        xml_path = require_file("bam_metadata_xml")
        summary = summarize_ena_xml(
            xml_path, label="ERZ6810918_BAM_metadata", output_dir=DIRS["validation"]
        )
        MANIFEST.log(ManifestEntry(
            step="parse_bam_metadata",
            description="Inspect ENA XML metadata for ERZ6810918 (BRIDGE panel "
                         "variants w/ MorexV3 reference) to determine actual file "
                         "content and whether raw read alignments are available",
            inputs=[str(xml_path)],
            outputs=list(summary["outputs"].values()),
            parameters={
                "n_total_elements": summary["n_total_elements"],
                "n_file_like_references": summary["n_file_like_references"],
                "n_accession_references": summary["n_accession_references"],
            },
            status="exploratory",
            notes="ERZ accessions are typically ENA 'Analysis' submissions. "
                  "IMPORTANT: the study description references 'variants', "
                  "suggesting this may be a VCF submission rather than raw "
                  "BAM/CRAM alignments. Check file_references output for actual "
                  "file extensions before assuming split-read/discordant-read "
                  "validation (Step 10) is possible from this source.",
        ))

    elif args.step == "parse-microbiome-metadata":
        xml_path = require_file("microbiome_metadata_xml")
        summary = summarize_ena_xml(
            xml_path, label="PRJEB50061_microbiome_metadata", output_dir=DIRS["microbiome"]
        )
        MANIFEST.log(ManifestEntry(
            step="parse_microbiome_metadata",
            description="Inspect ENA XML metadata for PRJEB50061 (Escudero-Martinez "
                         "et al. 2022 rhizosphere microbiome project) to determine "
                         "actual sample/run content and data availability",
            inputs=[str(xml_path)],
            outputs=list(summary["outputs"].values()),
            parameters={
                "n_total_elements": summary["n_total_elements"],
                "n_file_like_references": summary["n_file_like_references"],
                "n_accession_references": summary["n_accession_references"],
            },
            status="exploratory",
            notes="PRJEB50061 is a PROJECT-level accession with component "
                  "sub-projects (PRJEB50019/20/21/50400). A project-level XML "
                  "typically does NOT contain an OTU abundance table directly -- "
                  "that is normally a processed supplementary file from the "
                  "associated publication, or must be generated by running the "
                  "raw 16S/ITS amplicon reads through a pipeline (e.g. QIIME2/"
                  "DADA2). Check accession_references output for the actual "
                  "run/sample accessions referenced, to determine what further "
                  "retrieval (raw reads vs. processed table) is needed.",
    ))
    
    log.info("Step complete.")


if __name__ == "__main__":
    main()