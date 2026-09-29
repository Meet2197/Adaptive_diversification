
import pandas as pd
from pathlib import Path

BASE_DIR = Path(r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis")
ANNOT_PATH = BASE_DIR / "11_Functional_Annotation" / "SV_functional_annotation.csv"
CATALOGUE_PATH = BASE_DIR / "07_SV_Catalogue" / "nonredundant_SV_catalogue.csv"


def main():
    if not ANNOT_PATH.exists():
        print(f"NOT FOUND: {ANNOT_PATH} -- annotate step has not completed successfully yet.")
        return

    annot = pd.read_csv(ANNOT_PATH)
    catalogue = pd.read_csv(CATALOGUE_PATH)

    print("=" * 80)
    print("FUNCTIONAL ANNOTATION VERIFICATION")
    print("=" * 80)
    print(f"Total annotated SV records: {len(annot):,}")
    print(f"Total catalogue loci: {len(catalogue):,}")
    print(f"Match: {'YES' if len(annot) == len(catalogue) else 'NO -- MISMATCH, investigate'}")

    print(f"\n--- evidence_type breakdown ---")
    print(annot["evidence_type"].value_counts().to_string())

    n_gene_overlap = (annot["gene_overlap"].fillna("") != "").sum()
    n_exon_overlap = (annot["exon_overlap"].fillna("") != "").sum()
    n_promoter = (annot["promoter_proximity_genes"].fillna("") != "").sum()
    n_any_hit = (
        (annot["gene_overlap"].fillna("") != "") |
        (annot["exon_overlap"].fillna("") != "") |
        (annot["promoter_proximity_genes"].fillna("") != "")
    ).sum()

    print(f"\n--- Overlap statistics ---")
    print(f"SVs with gene-body overlap:       {n_gene_overlap:,} ({100*n_gene_overlap/len(annot):.2f}%)")
    print(f"SVs with exon overlap:            {n_exon_overlap:,} ({100*n_exon_overlap/len(annot):.2f}%)")
    print(f"SVs with promoter proximity hit:  {n_promoter:,} ({100*n_promoter/len(annot):.2f}%)")
    print(f"SVs with ANY functional evidence: {n_any_hit:,} ({100*n_any_hit/len(annot):.2f}%)")

    if n_any_hit == 0:
        print("\n*** WARNING: ZERO overlaps found across all SVs. This strongly ***")
        print("*** suggests a remaining chromosome-naming or coordinate issue. ***")
        print("*** Do NOT trust this output -- investigate before proceeding.  ***")
    else:
        print(f"\nResult looks genuine -- {n_any_hit:,} SVs show real positional evidence.")

    # Sample a few real hits for manual eyeballing
    hits = annot[annot["gene_overlap"].fillna("") != ""]
    if len(hits) > 0:
        print(f"\n--- Sample gene-overlap hits (first 5) ---")
        print(hits[["SV_ID", "gene_overlap", "exon_overlap"]].head(5).to_string(index=False))

    annot.to_csv(BASE_DIR / "11_Functional_Annotation" / "SV_functional_annotation_VERIFIED.csv", index=False)


if __name__ == "__main__":
    main()