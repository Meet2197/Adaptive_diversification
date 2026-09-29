
from pathlib import Path
import gzip
import re
from collections import Counter

GFF3_PATH = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis\external_data\genes.gff3")

# Your SV catalogue's confirmed chromosome naming (from earlier validated QC)
CATALOGUE_CHROMOSOMES = {f"chr{i}H" for i in range(1, 8)} | {"chrUn"}


def open_maybe_gzip(path: Path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def inspect():
    if not GFF3_PATH.exists():
        raise FileNotFoundError(f"{GFF3_PATH} not found.")

    print("=" * 80)
    print("GFF3 FILE INSPECTION")
    print("=" * 80)
    print(f"File: {GFF3_PATH}")
    print(f"Size: {GFF3_PATH.stat().st_size:,} bytes\n")

    header_lines = []
    seqid_counter = Counter()
    featuretype_counter = Counter()
    n_data_lines = 0
    n_malformed = 0
    min_coords_by_seqid = {}
    max_coords_by_seqid = {}
    sample_gene_lines = []

    with open_maybe_gzip(GFF3_PATH) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#"):
                if len(header_lines) < 15:
                    header_lines.append(line.rstrip())
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                n_malformed += 1
                continue

            seqid, source, featuretype, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
            n_data_lines += 1
            seqid_counter[seqid] += 1
            featuretype_counter[featuretype] += 1

            try:
                start_i, end_i = int(start), int(end)
                if seqid not in min_coords_by_seqid or start_i < min_coords_by_seqid[seqid]:
                    min_coords_by_seqid[seqid] = start_i
                if seqid not in max_coords_by_seqid or end_i > max_coords_by_seqid[seqid]:
                    max_coords_by_seqid[seqid] = end_i
            except ValueError:
                n_malformed += 1

            if featuretype == "gene" and len(sample_gene_lines) < 3:
                sample_gene_lines.append(line.rstrip())

            if line_num % 500_000 == 0:
                print(f"  ...processed {line_num:,} lines so far")

    print("\n--- Header / pragma lines (first 15) ---")
    for h in header_lines:
        print(f"  {h}")

    print(f"\n--- Summary ---")
    print(f"Total data lines: {n_data_lines:,}")
    print(f"Malformed lines: {n_malformed:,}")

    print(f"\n--- Unique seqid (chromosome) names found in GFF3 ({len(seqid_counter)}) ---")
    for seqid, count in seqid_counter.most_common(20):
        span = f"{min_coords_by_seqid.get(seqid, '?'):,} - {max_coords_by_seqid.get(seqid, '?'):,}"
        print(f"  {seqid:20s}  {count:>10,} features   range: {span}")
    if len(seqid_counter) > 20:
        print(f"  ... and {len(seqid_counter) - 20} more")

    print(f"\n--- Unique feature types (column 3) ---")
    for ftype, count in featuretype_counter.most_common():
        print(f"  {ftype:20s}  {count:>10,}")

    print(f"\n--- Sample 'gene' feature lines (raw) ---")
    for line in sample_gene_lines:
        print(f"  {line}")

    # ------------------------------------------------------------
    # CRITICAL CHECK: does GFF3 seqid naming match the SV catalogue?
    # ------------------------------------------------------------
    print("\n" + "=" * 80)
    print("CHROMOSOME NAMING COMPATIBILITY CHECK (critical)")
    print("=" * 80)

    gff3_seqids = set(seqid_counter.keys())
    exact_match = gff3_seqids & CATALOGUE_CHROMOSOMES
    print(f"GFF3 seqids: {sorted(gff3_seqids)[:15]}{'...' if len(gff3_seqids) > 15 else ''}")
    print(f"Catalogue chromosomes: {sorted(CATALOGUE_CHROMOSOMES)}")
    print(f"\nExact-name matches: {len(exact_match)} / {len(CATALOGUE_CHROMOSOMES)}")

    if len(exact_match) == len(CATALOGUE_CHROMOSOMES):
        print("\n=> RESULT: GFF3 seqid naming MATCHES the SV catalogue exactly. "
              "No harmonization needed -- safe to proceed directly to "
              "gffutils database creation and annotate_svs().")
    else:
        print("\n=> RESULT: MISMATCH DETECTED. GFF3 seqid naming does NOT match "
              "the SV catalogue's 'chr1H'..'chr7H'/'chrUn' convention. "
              "A chromosome-name mapping/harmonization step is REQUIRED "
              "before running annotate_svs(), otherwise every gene-overlap "
              "lookup will silently return empty results.")

        # Try to suggest a mapping heuristic
        print("\nAttempting automatic mapping suggestion:")
        for cat_chrom in sorted(CATALOGUE_CHROMOSOMES):
            core = cat_chrom.replace("chr", "")  # e.g. "1H", "Un"
            candidates = [s for s in gff3_seqids if core.lower() in s.lower()]
            print(f"  {cat_chrom:8s} -> candidate GFF3 seqid(s): {candidates[:3]}")

    # Feature-type check
    print("\n" + "=" * 80)
    print("FEATURE TYPE COMPATIBILITY CHECK")
    print("=" * 80)
    for required in ["gene", "exon"]:
        present = required in featuretype_counter
        print(f"  '{required}': {'FOUND' if present else 'NOT FOUND'} "
              f"({featuretype_counter.get(required, 0):,} features)")


if __name__ == "__main__":
    inspect()