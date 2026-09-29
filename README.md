

# Barley Structural Variant Catalogue and Pangenomics Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.11-blue.svg)

A reproducible pipeline for constructing a non-redundant structural
variant (SV) catalogue across a diversity panel of barley (*Hordeum
vulgare* L.) accessions from Assemblytics assembly-comparison output,
with single-linkage chain-artifact correction, functional annotation
against the MorexV3 gene model, and accession-provenance resolution.

---

## What this repository does

Starting from 74 whole-genome assembly comparisons against the MorexV3
reference (Assemblytics output), this pipeline:

1. Validates and cleans raw structural variant (SV) calls
2. Normalizes SV type nomenclature and documents coordinate semantics
   distinct between breakpoint-type and repeat-associated SVs
3. Clusters equivalent SVs across accessions using reciprocal-overlap
   and size-ratio criteria, with a purpose-built correction for
   single-linkage clustering artifacts (see [Methodological Note](#methodological-note))
4. Builds a non-redundant, frequency-classified SV catalogue
5. Constructs an accession × SV genotype matrix
6. Annotates candidate SVs against the MorexV3 IPK gene model (gene
   body, exon, and promoter-proximity overlap)
7. Resolves accession identity against IPK genebank and BRIDGE panel
   passport records

---

## Repository structure

## Installation

```bash
conda create -n barley_sv python=3.11
conda activate barley_sv
pip install -r requirements.txt
```

`requirements.txt`:

<pre node="[object Object]"><div class="relative block bg-neutral-800 rounded-lg p-1 my-4"></div></pre>

```Python
pandas
numpy
scipy
matplotlib
seaborn
gffutils
statsmodels
requests
openpyxl
```

## Usage

```Python
# Step 1: parse, clean, and QC raw Assemblytics BED files
python main.py

# Step 2-3: scalable clustering with chain-split correction
python extended_analysis.py cluster --max-distance 500 --min-overlap 0.7

# Step 4: build the non-redundant SV catalogue
python extended_analysis.py catalogue

# Step 5: build the accession x SV genotype matrix
python extended_analysis.py genotype-matrix

# Step 6: functional annotation against MorexV3 gene model
python extended_analysis.py annotate
python verify_annotation_results.py

# Accession crosswalk
python build_accession_crosswalk.py
python crosswalk_progress_check.py
```

Input data (Assemblytics `.bed` files, MorexV3 GFF3 annotation) must be placed in the paths configured at the top of `main.py` / `extended_analysis.py`, or the scripts will raise an explicit `FileNotFoundError` rather than proceeding silently.

## Key results

| Metric                                              | Value           |
| --------------------------------------------------- | --------------- |
| Accessions analyzed                                 | 74              |
| Raw cleaned SV records                              | 8,396,644       |
| Non-redundant SV loci (post chain-split correction) | 1,703,270       |
| Loci with gene-body overlap                         | 122,662 (7.20%) |
| Loci with exon overlap                              | 97,582 (5.73%)  |
| Loci with any functional evidence                   | 167,120 (9.81%) |
| Genes in annotation reference                       | 34,706          |
| Accessions with confirmed genebank identity         | 36 / 74 (48.6%) |

## Methodological note

A single-linkage connected-components clustering step, when applied naively at genome scale, is vulnerable to a "chaining" artifact: SV records that individually satisfy pairwise merge criteria only with immediate neighbors can become transitively linked into clusters spanning megabases, despite no single pair of members within that cluster actually meeting the merge threshold directly. This was identified empirically (initial genome-wide run produced clusters with genomic spans up to 3.19 Mb) and corrected using a post-hoc, diameter-constrained re-splitting step, validated identically at both chromosome-scale (930,598 records) and genome-scale (8,396,644 records). See `extended_analysis.py::split_chained_clusters()` for the implementation and inline documentation.

---

## Data availability

Raw Assemblytics output and the MorexV3 assembly/annotation are distributed by IPK Gatersleben via the e!DAL repository. Accession identity cross-referencing used passport data from Milner et al. (2019),  *Nature Genetics* , DOI: [10.5447/IPK/2018/9](https://doi.org/10.5447/IPK/2018/9).

## Citation

If you use this pipeline, please cite this repository and the associated data sources listed above.
