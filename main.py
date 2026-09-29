from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import re
import warnings

warnings.filterwarnings("ignore")

# ============================================================
# CONFIGURATION
# ============================================================

INPUT_DIR = Path(
    r"C:\Users\bhatt\Desktop\Omics\diversification_barley\Assemblytics"
)

OUTPUT_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")
MIN_SV_SIZE = 50
MAX_MAIN_SV_SIZE = 1_000_000
SV_CLUSTER_BIN = 1000
MIN_ACCESSION_SUPPORT = 2
KNOWN_PLACED_CHROMOSOMES = {f"chr{i}H" for i in range(1, 8)}

QC_DIR = OUTPUT_DIR / "01_QC"
SUMMARY_DIR = OUTPUT_DIR / "02_SV_Summary"
SHARED_DIR = OUTPUT_DIR / "03_Shared_Private_SVs"
MATRIX_DIR = OUTPUT_DIR / "04_PanGWAS_Matrix"
FIGURE_DIR = OUTPUT_DIR / "05_Figures"
FINAL_DIR = OUTPUT_DIR / "06_Final_Datasets"

for directory in [
    QC_DIR,
    SUMMARY_DIR,
    SHARED_DIR,
    MATRIX_DIR,
    FIGURE_DIR,
    FINAL_DIR
]:
    directory.mkdir(
        parents=True,
        exist_ok=True
    )


# ============================================================
# ACCESSION NAME
# ============================================================

def accession_from_filename(path):

    suffix = ".Assemblytics_structural_variants.bed"

    if path.name.endswith(suffix):

        return path.name[:-len(suffix)]

    return path.stem


# ============================================================
# ASSEMBLYTICS READER
# ============================================================

def read_assemblytics(path):

    accession = accession_from_filename(path)

    expected_columns = [
        "reference",
        "ref_start",
        "ref_stop",
        "ID",
        "size",
        "strand",
        "type",
        "ref_gap_size",
        "query_gap_size",
        "query_coordinates",
        "method"
    ]

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=expected_columns,
        dtype=str,
        encoding="utf-8-sig"
    )

    if df.empty:

        raise ValueError(
            f"{path.name}: no SV records found"
        )

    # --------------------------------------------------------
    # Convert numeric columns
    # --------------------------------------------------------

    numeric_columns = [
        "ref_start",
        "ref_stop",
        "size",
        "ref_gap_size",
        "query_gap_size"
    ]

    for col in numeric_columns:

        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )

    # --------------------------------------------------------
    # Remove invalid records
    # --------------------------------------------------------

    original_records = len(df)

    invalid_mask = (
        df["reference"].isna()
        |
        df["ref_start"].isna()
        |
        df["ref_stop"].isna()
        |
        df["size"].isna()
        |
        df["type"].isna()
    )

    invalid_records = invalid_mask.sum()

    df = df.loc[
        ~invalid_mask
    ].copy()

    # --------------------------------------------------------
    # Coordinate validation
    # --------------------------------------------------------

    bad_coordinates = (
        (df["ref_start"] < 0)
        |
        (df["ref_stop"] < df["ref_start"])
    )

    invalid_coordinate_count = (
        bad_coordinates.sum()
    )

    df = df.loc[
        ~bad_coordinates
    ].copy()

    # --------------------------------------------------------
    # Size validation
    # --------------------------------------------------------

    df["SV_SIZE"] = df["size"].abs()

    df = df[
        df["SV_SIZE"] >= MIN_SV_SIZE
    ].copy()

    # --------------------------------------------------------
    # Normalize SV type
    # --------------------------------------------------------

    df["SV_TYPE"] = (
        df["type"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    type_map = {

        "insertion": "Insertion",
        "insertions": "Insertion",

        "deletion": "Deletion",
        "deletions": "Deletion",

        "tandem_expansion":
            "Tandem_expansion",

        "tandem_contraction":
            "Tandem_contraction",

        "repeat_expansion":
            "Repeat_expansion",

        "repeat_contraction":
            "Repeat_contraction"
    }

    df["SV_TYPE_NORMALIZED"] = (
        df["SV_TYPE"]
        .map(type_map)
        .fillna(
            df["SV_TYPE"].str.title()
        )
    )

    df["REF_SPAN_BP"] = df["ref_stop"] - df["ref_start"]
    df["EXCEEDS_MAX_PLAUSIBLE_SIZE"] = (
        (df["SV_SIZE"] > MAX_MAIN_SV_SIZE) |
        (df["REF_SPAN_BP"] > MAX_MAIN_SV_SIZE)
    )

    df["accession"] = accession

    df["source_file"] = str(
        path.resolve()
    )

    df["SV_RECORD_ID"] = [
        f"{accession}_SV_{i+1:08d}"
        for i in range(len(df))
    ]

    qc = {
        "accession": accession,
        "file": str(path),
        "original_records": original_records,
        "invalid_records": int(
            invalid_records
        ),
        "invalid_coordinates": int(
            invalid_coordinate_count
        ),
        "final_records": len(df),
        "chromosomes": df[
            "reference"
        ].nunique(),
        "min_size": df[
            "SV_SIZE"
        ].min(),
        "median_size": df[
            "SV_SIZE"
        ].median(),
        "mean_size": df[
            "SV_SIZE"
        ].mean(),
        "max_size": df[
            "SV_SIZE"
        ].max(),
        "oversized_flagged": int(
            df["EXCEEDS_MAX_PLAUSIBLE_SIZE"].sum()
        )
    }

    return df, qc


# ============================================================
# STEP 1
# DATASET DISCOVERY + QC
# ============================================================

print("\n")
print("=" * 80)
print("STEP 1 — DATASET DISCOVERY AND QUALITY CONTROL")
print("=" * 80)

files = sorted(
    INPUT_DIR.rglob(
        "*Assemblytics_structural_variants.bed"
    )
)

print(
    f"\nAssemblytics files found: {len(files)}"
)

if len(files) == 0:

    raise SystemExit(
        "\nERROR: No Assemblytics BED files found."
    )


datasets = []
qc_records = []
failed_files = []

for file in files:

    print(
        f"Reading: {file.name}"
    )

    try:

        df, qc = read_assemblytics(
            file
        )

        datasets.append(df)
        qc_records.append(qc)

        print(
            f"  OK: {len(df):,} SVs"
        )

    except Exception as e:

        print(
            f"  ERROR: {e}"
        )

        failed_files.append({
            "file": str(file),
            "error": str(e)
        })


# ------------------------------------------------------------
# Stop safely if no files succeeded
# ------------------------------------------------------------

if not datasets:

    print("\n")
    print("=" * 80)
    print("FATAL ERROR")
    print("=" * 80)
    print(
        "No Assemblytics files could be parsed."
    )

    pd.DataFrame(
        failed_files
    ).to_csv(
        QC_DIR / "failed_files.csv",
        index=False
    )

    raise SystemExit(1)


# ------------------------------------------------------------
# Save QC
# ------------------------------------------------------------

qc_df = pd.DataFrame(
    qc_records
)

qc_df.to_csv(
    QC_DIR / "dataset_QC_report.csv",
    index=False
)

if failed_files:

    pd.DataFrame(
        failed_files
    ).to_csv(
        QC_DIR / "failed_files.csv",
        index=False
    )


# ============================================================
# COMBINE DATASETS
# ============================================================

all_sv = pd.concat(
    datasets,
    ignore_index=True
)

print(
    f"\nTotal cleaned SV records: "
    f"{len(all_sv):,}"
)

print(
    f"Accessions: "
    f"{all_sv['accession'].nunique():,}"
)

print(
    f"Reference chromosomes: "
    f"{all_sv['reference'].nunique():,}"
)

chrom_composition = all_sv["reference"].value_counts().reset_index()
chrom_composition.columns = ["reference", "record_count"]
chrom_composition["is_placed_chromosome"] = chrom_composition["reference"].isin(
    KNOWN_PLACED_CHROMOSOMES
)
chrom_composition.to_csv(QC_DIR / "chromosome_composition_report.csv", index=False)

print("\nChromosome composition:")
print(chrom_composition.to_string(index=False))

unplaced_count = chrom_composition.loc[
    ~chrom_composition["is_placed_chromosome"], "record_count"
].sum()
if unplaced_count > 0:
    pct = 100 * unplaced_count / len(all_sv)
    print(
        f"\nNOTE: {unplaced_count:,} records ({pct:.2f}%) mapped to unplaced "
        f"scaffold sequence. Retained in dataset; excluded from per-chromosome "
        f"figures and positional/functional annotation requiring placed coordinates."
    )

oversized_report = all_sv[all_sv["EXCEEDS_MAX_PLAUSIBLE_SIZE"]][
    ["accession", "reference", "ref_start", "ref_stop", "REF_SPAN_BP", "SV_SIZE",
     "SV_TYPE_NORMALIZED", "SV_RECORD_ID"]
]

oversized_report.to_csv(QC_DIR / "oversized_SV_records_flagged.csv", index=False)
print(
    f"\nFlagged {len(oversized_report):,} SV records exceeding "
    f"{MAX_MAIN_SV_SIZE:,}bp (retained in dataset, listed in "
    f"oversized_SV_records_flagged.csv for manual review)"
)

all_sv["SIZE_CLASS"] = pd.cut(

    all_sv["SV_SIZE"],

    bins=[
        0,
        100,
        500,
        1000,
        5000,
        10000,
        50000,
        100000,
        np.inf
    ],

    labels=[
        "50–100 bp",
        "101–500 bp",
        "501 bp–1 kb",
        "1–5 kb",
        "5–10 kb",
        "10–50 kb",
        "50–100 kb",
        ">100 kb"
    ],

    include_lowest=True
)


# ============================================================
# SAVE CLEANED MASTER DATASET
# ============================================================

all_sv.to_csv(
    FINAL_DIR /
    "ALL_Assemblytics_SVs_CLEANED.csv",
    index=False
)


# ============================================================
# STEP 2
# SV QUANTIFICATION
# ============================================================

print("\n")
print("=" * 80)
print("STEP 2 — SV QUANTIFICATION AND CHARACTERIZATION")
print("=" * 80)


# ------------------------------------------------------------
# Accession summary
# ------------------------------------------------------------

accession_summary = (

    all_sv

    .groupby("accession")

    .agg(

        total_SVs=(
            "SV_RECORD_ID",
            "count"
        ),

        total_SV_bp=(
            "SV_SIZE",
            "sum"
        ),

        mean_SV_size=(
            "SV_SIZE",
            "mean"
        ),

        median_SV_size=(
            "SV_SIZE",
            "median"
        ),

        largest_SV=(
            "SV_SIZE",
            "max"
        ),

        chromosomes=(
            "reference",
            "nunique"
        )
    )

    .reset_index()
)

accession_summary.to_csv(
    SUMMARY_DIR /
    "accession_SV_summary.csv",
    index=False
)


# ------------------------------------------------------------
# SV type summary
# ------------------------------------------------------------

sv_type_summary = (

    all_sv

    .groupby(
        [
            "accession",
            "SV_TYPE_NORMALIZED"
        ]
    )

    .size()

    .reset_index(
        name="SV_count"
    )
)

sv_type_summary.to_csv(
    SUMMARY_DIR /
    "SV_type_by_accession.csv",
    index=False
)


# ------------------------------------------------------------
# Size distribution
# ------------------------------------------------------------

size_summary = (

    all_sv

    .groupby(
        "SIZE_CLASS",
        observed=False
    )

    .size()

    .reset_index(
        name="SV_count"
    )
)

size_summary.to_csv(
    SUMMARY_DIR /
    "SV_size_distribution.csv",
    index=False
)


# ------------------------------------------------------------
# Chromosome distribution
# ------------------------------------------------------------

chromosome_summary = (

    all_sv

    .groupby(
        [
            "accession",
            "reference"
        ]
    )

    .size()

    .reset_index(
        name="SV_count"
    )
)

chromosome_summary.to_csv(
    SUMMARY_DIR /
    "SV_by_chromosome.csv",
    index=False
)


# ------------------------------------------------------------
# Insertion/deletion summary
# ------------------------------------------------------------

indel_summary = (

    all_sv

    .groupby(
        [
            "accession",
            "SV_TYPE_NORMALIZED"
        ]
    )

    .agg(
        count=(
            "SV_RECORD_ID",
            "count"
        ),

        total_bp=(
            "SV_SIZE",
            "sum"
        ),

        median_size=(
            "SV_SIZE",
            "median"
        )
    )

    .reset_index()
)

indel_summary.to_csv(
    SUMMARY_DIR /
    "SV_type_size_summary.csv",
    index=False
)


# ============================================================
# STEP 3
# SHARED / PRIVATE SV ANALYSIS
# ============================================================

print("\n")
print("=" * 80)
print("STEP 3 — SHARED AND PRIVATE SV ANALYSIS")
print("=" * 80)


all_sv["COORDINATE_BIN"] = (

    all_sv["ref_start"]
    // SV_CLUSTER_BIN
).astype(int)


all_sv["SV_SIGNATURE"] = (

    all_sv["reference"]
    .astype(str)

    + "|"

    + all_sv["SV_TYPE_NORMALIZED"]
    .astype(str)

    + "|"

    + all_sv["COORDINATE_BIN"]
    .astype(str)
)


# ------------------------------------------------------------
# Support
# ------------------------------------------------------------

signature_support = (

    all_sv

    .groupby(
        "SV_SIGNATURE"
    )

    .agg(

        accession_count=(
            "accession",
            "nunique"
        ),

        record_count=(
            "SV_RECORD_ID",
            "count"
        ),

        median_size=(
            "SV_SIZE",
            "median"
        ),

        reference=(
            "reference",
            "first"
        ),

        SV_type=(
            "SV_TYPE_NORMALIZED",
            "first"
        )
    )

    .reset_index()
)


signature_support.to_csv(
    SHARED_DIR /
    "SV_signature_support.csv",
    index=False
)


# ------------------------------------------------------------
# Classify
# ------------------------------------------------------------

all_sv = all_sv.merge(

    signature_support[
        [
            "SV_SIGNATURE",
            "accession_count"
        ]
    ],

    on="SV_SIGNATURE",

    how="left"
)


all_sv["SHARING_CLASS"] = np.where(

    all_sv["accession_count"] >= 2,

    "Shared",

    "Private"
)


shared_sv = all_sv[
    all_sv["SHARING_CLASS"] == "Shared"
].copy()


private_sv = all_sv[
    all_sv["SHARING_CLASS"] == "Private"
].copy()


shared_sv.to_csv(
    SHARED_DIR /
    "shared_SVs_exploratory.csv",
    index=False
)


private_sv.to_csv(
    SHARED_DIR /
    "private_SVs_exploratory.csv",
    index=False
)


# ============================================================
# STEP 4
# PanGWAS MATRIX (EXPLORATORY — see end-of-script warning)
# ============================================================

print("\n")
print("=" * 80)
print("STEP 4 — EXPLORATORY SV PRESENCE/ABSENCE MATRIX")
print("=" * 80)

sv_matrix = pd.crosstab(

    all_sv[
        "SV_SIGNATURE"
    ],

    all_sv[
        "accession"
    ]

)

sv_matrix = (

    sv_matrix
    .clip(upper=1)
    .astype(int)
)


sv_matrix.to_csv(
    MATRIX_DIR /
    "SV_presence_absence_ALL_EXPLORATORY.csv"
)


# ------------------------------------------------------------
# High-confidence shared matrix (still exploratory)
# ------------------------------------------------------------

high_confidence = (

    signature_support[
        signature_support[
            "accession_count"
        ]
        >= MIN_ACCESSION_SUPPORT
    ]

    [
        "SV_SIGNATURE"
    ]
)


sv_matrix_high = sv_matrix.loc[
    sv_matrix.index.isin(
        high_confidence
    )
]


sv_matrix_high.to_csv(
    MATRIX_DIR /
    "SV_presence_absence_HIGH_CONFIDENCE_EXPLORATORY.csv"
)


# ------------------------------------------------------------
# SV metadata
# ------------------------------------------------------------

sv_metadata = (

    all_sv

    .groupby(
        "SV_SIGNATURE"
    )

    .agg(

        reference=(
            "reference",
            "first"
        ),

        start=(
            "ref_start",
            "min"
        ),

        end=(
            "ref_stop",
            "max"
        ),

        SV_type=(
            "SV_TYPE_NORMALIZED",
            "first"
        ),

        median_size=(
            "SV_SIZE",
            "median"
        ),

        min_size=(
            "SV_SIZE",
            "min"
        ),

        max_size=(
            "SV_SIZE",
            "max"
        ),

        accession_count=(
            "accession_count",
            "first"
        )
    )

    .reset_index()
)


sv_metadata.to_csv(
    MATRIX_DIR /
    "SV_metadata_ALL_EXPLORATORY.csv",
    index=False
)


# ------------------------------------------------------------
# Exploratory PanGWAS-ready package (DO NOT use for final analysis —
# see end-of-script warning)
# ------------------------------------------------------------

PAN_DIR = FINAL_DIR / "PanGWAS_EXPLORATORY_DO_NOT_USE_FOR_FINAL"

PAN_DIR.mkdir(
    parents=True,
    exist_ok=True
)


sv_matrix_high.to_csv(
    PAN_DIR /
    "barley_SV_PanGWAS_presence_absence_EXPLORATORY.csv"
)


sv_metadata[
    sv_metadata[
        "SV_SIGNATURE"
    ].isin(
        sv_matrix_high.index
    )
].to_csv(
    PAN_DIR /
    "barley_SV_PanGWAS_metadata_EXPLORATORY.csv",
    index=False
)


accession_summary.to_csv(
    PAN_DIR /
    "barley_accession_summary.csv",
    index=False
)


# ============================================================
# STEP 5
# FIGURES
# ============================================================

print("\n")
print("=" * 80)
print("STEP 5 — PUBLICATION-READY FIGURES")
print("=" * 80)

plot_df = accession_summary.sort_values(
    "total_SVs",
    ascending=False
)


plt.figure(
    figsize=(16, 7)
)

sns.barplot(
    data=plot_df,
    x="accession",
    y="total_SVs"
)

plt.xticks(
    rotation=90,
    fontsize=7
)

plt.xlabel(
    "Barley accession"
)

plt.ylabel(
    "Number of SVs"
)

plt.title(
    "Assemblytics Structural Variant Counts Across Barley Accessions"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR /
    "Figure1_SV_count_by_accession.png",
    dpi=400
)

plt.close()


# ------------------------------------------------------------
# Figure 2
# SV types
# ------------------------------------------------------------

type_pivot = pd.crosstab(
    all_sv["accession"],
    all_sv["SV_TYPE_NORMALIZED"]
)


type_pivot.plot(
    kind="bar",
    stacked=True,
    figsize=(16, 7)
)

plt.xticks(
    rotation=90,
    fontsize=7
)

plt.xlabel(
    "Barley accession"
)

plt.ylabel(
    "SV count"
)

plt.title(
    "Structural Variant Classes Across Barley Accessions"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR /
    "Figure2_SV_type_distribution.png",
    dpi=400
)

plt.close()


# ------------------------------------------------------------
# Figure 3
# SV size
# ------------------------------------------------------------

plt.figure(
    figsize=(12, 7)
)

sns.histplot(
    data=all_sv,
    x="SV_SIZE",
    hue="SV_TYPE_NORMALIZED",
    bins=100,
    element="step",
    log_scale=True
)

plt.xlabel(
    "SV size (bp)"
)

plt.ylabel(
    "Number of SVs"
)

plt.title(
    "Assemblytics Structural Variant Size Distribution"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR /
    "Figure3_SV_size_distribution.png",
    dpi=400
)

plt.close()


# ------------------------------------------------------------
# Figure 4
# Chromosome distribution (placed chromosomes only)
# ------------------------------------------------------------

placed_sv = all_sv[all_sv["reference"].isin(KNOWN_PLACED_CHROMOSOMES)]

chr_pivot = pd.crosstab(
    placed_sv["reference"],
    placed_sv["accession"]
)


plt.figure(
    figsize=(16, 8)
)

sns.heatmap(
    chr_pivot,
    cmap="viridis"
)

plt.xlabel(
    "Barley accession"
)

plt.ylabel(
    "MorexV3 reference chromosome"
)

plt.title(
    "Chromosome-Level Distribution of Structural Variants\n(placed chromosomes only; chrUn excluded)"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR /
    "Figure4_chromosome_distribution.png",
    dpi=400
)

plt.close()


# ------------------------------------------------------------
# Figure 5
# SV sharing (exploratory signature basis)
# ------------------------------------------------------------

sharing = (

    signature_support

    .groupby(
        "accession_count"
    )

    .size()

    .reset_index(
        name="SV_count"
    )
)


plt.figure(
    figsize=(11, 7)
)

sns.barplot(
    data=sharing,
    x="accession_count",
    y="SV_count"
)

plt.xlabel(
    "Number of accessions containing SV"
)

plt.ylabel(
    "Number of consolidated SV signatures (exploratory)"
)

plt.title(
    "Structural Variant Sharing Across Barley Accessions (Exploratory Signature)"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR /
    "Figure5_SV_sharing_EXPLORATORY.png",
    dpi=400
)

plt.close()


# ------------------------------------------------------------
# Figure 6
# SV-based accession clustering (exploratory signature basis)
# ------------------------------------------------------------

if sv_matrix_high.shape[0] >= 2:

    clustering_matrix = (
        sv_matrix_high.T
    )

    if clustering_matrix.shape[0] >= 3:

        distance = pdist(
            clustering_matrix,
            metric="jaccard"
        )

        linkage_matrix = linkage(
            distance,
            method="average"
        )

        plt.figure(
            figsize=(16, 8)
        )

        dendrogram(
            linkage_matrix,
            labels=clustering_matrix.index,
            leaf_rotation=90
        )

        plt.xlabel(
            "Barley accession"
        )

        plt.ylabel(
            "Jaccard distance"
        )

        plt.title(
            "Barley Accession Clustering Based on SV Presence/Absence (Exploratory)"
        )

        plt.tight_layout()

        plt.savefig(
            FIGURE_DIR /
            "Figure6_SV_Jaccard_clustering_EXPLORATORY.png",
            dpi=400
        )

        plt.close()


# ============================================================
# MASTER SUMMARY
# ============================================================

master_summary = pd.DataFrame({

    "Metric": [

        "Assemblytics files discovered",

        "Successfully parsed files",

        "Failed files",

        "Barley accessions",

        "Total cleaned SV records",

        "Reference chromosomes (incl. chrUn)",

        "Placed chromosomes",

        "Unique SV signatures (exploratory)",

        "Shared SV records (exploratory)",

        "Private SV records (exploratory)",

        "High-confidence SV signatures (exploratory)",

        "Oversized SVs flagged (>MAX_MAIN_SV_SIZE)",

        "Minimum SV size",

        "Maximum SV size"

    ],

    "Value": [

        len(files),

        len(datasets),

        len(failed_files),

        all_sv[
            "accession"
        ].nunique(),

        len(all_sv),

        all_sv[
            "reference"
        ].nunique(),

        len(KNOWN_PLACED_CHROMOSOMES),

        all_sv[
            "SV_SIGNATURE"
        ].nunique(),

        len(shared_sv),

        len(private_sv),

        len(sv_matrix_high),

        int(all_sv["EXCEEDS_MAX_PLAUSIBLE_SIZE"].sum()),

        all_sv[
            "SV_SIZE"
        ].min(),

        all_sv[
            "SV_SIZE"
        ].max()
    ]
})


master_summary.to_csv(
    OUTPUT_DIR /
    "MASTER_Prompt1_Barley_SV_Summary.csv",
    index=False
)


# ============================================================
# GRAPH-PANGENOME READINESS REPORT
# ============================================================

report = f"""

Input:
{INPUT_DIR}

Assemblytics files discovered:
{len(files)}

Successfully parsed:
{len(datasets)}

Failed:
{len(failed_files)}

Barley accessions:
{all_sv["accession"].nunique()}

Total cleaned Assemblytics SV records:
{len(all_sv)}

Reference chromosomes (incl. chrUn):
{all_sv["reference"].nunique()}

Placed chromosomes:
{len(KNOWN_PLACED_CHROMOSOMES)}

Unique EXPLORATORY SV signatures:
{all_sv["SV_SIGNATURE"].nunique()}

Exploratory shared SV records:
{len(shared_sv)}

Exploratory private SV records:
{len(private_sv)}

High-confidence SV signatures (exploratory):
{len(sv_matrix_high)}

Oversized SVs flagged (> {MAX_MAIN_SV_SIZE:,} bp):
{int(all_sv["EXCEEDS_MAX_PLAUSIBLE_SIZE"].sum())}

Minimum SV size:
{all_sv["SV_SIZE"].min()} bp

Maximum SV size:
{all_sv["SV_SIZE"].max()} bp

"""

with open(
    OUTPUT_DIR /
    "GRAPH_PANGENOME_READINESS_REPORT.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(report)


# ============================================================
# FINISHED
# ============================================================

print("\n")
print("=" * 80)
print("PROMPT 1 — BARLEY ANALYSIS COMPLETED")
print("=" * 80)

print(
    f"""
Files discovered:              {len(files)}
Files successfully parsed:     {len(datasets)}
Failed files:                  {len(failed_files)}
Barley accessions:             {all_sv['accession'].nunique():,}
Cleaned SV records:            {len(all_sv):,}
SV signatures (exploratory):   {all_sv['SV_SIGNATURE'].nunique():,}
High-confidence (exploratory): {len(sv_matrix_high):,}
Oversized SVs flagged:         {int(all_sv['EXCEEDS_MAX_PLAUSIBLE_SIZE'].sum()):,}

OUTPUT:
{OUTPUT_DIR}

Exploratory PanGWAS matrix:
{PAN_DIR / 'barley_SV_PanGWAS_presence_absence_EXPLORATORY.csv'}

Exploratory SV metadata:
{PAN_DIR / 'barley_SV_PanGWAS_metadata_EXPLORATORY.csv'}

Master summary:
{OUTPUT_DIR / 'MASTER_Prompt1_Barley_SV_Summary.csv'}

Readiness report:
{OUTPUT_DIR / 'GRAPH_PANGENOME_READINESS_REPORT.txt'}
"""
)

print("\n" + "!" * 80)
print("IMPORTANT: The PanGWAS matrix produced by THIS script is EXPLORATORY")
print("(built from a 1kb coordinate-bin heuristic, not reciprocal-overlap")
print("clustering). DO NOT use it for final SV-PanGWAS or the manuscript.")
print("Run extended_analysis.py (cluster -> catalogue -> genotype-matrix)")
print("to produce the validated, chain-split-corrected non-redundant SV")
print("catalogue and genotype matrix in 07_SV_Catalogue/ and")
print("08_Genotype_Matrix/ instead.")
print("!" * 80)