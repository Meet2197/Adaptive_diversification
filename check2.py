import pandas as pd
import numpy as np

df = pd.read_csv(
    r"C:\Users\bhatt\Desktop\Omics\Adaptive_diversification\Barley_Analysis\07_SV_Catalogue\clustered_SVs.csv"
)

summary = df.groupby("CLUSTER_ID").agg(
    n_records=("SV_RECORD_ID", "count"),
    min_start=("ref_start", "min"),
    max_start=("ref_start", "max"),      # NEW: needed to check TRUE chaining
    max_end=("ref_stop", "max"),
    any_chain_split=("CHAIN_SPLIT_APPLIED", "any"),
)

# The metric the chain-split fix ACTUALLY controls:
summary["start_diameter"] = summary["max_start"] - summary["min_start"]

# The (potentially misleading) metric used in the previous diagnostic:
summary["span_bp_old_metric"] = summary["max_end"] - summary["min_start"]

multi = summary[summary["n_records"] > 1].copy()

print("=== TRUE chaining check: start_diameter (what the fix actually controls) ===")
print(multi["start_diameter"].describe())
print(f"\nMulti-record clusters with start_diameter > 500bp (should be ~0): "
      f"{(multi['start_diameter'] > 500).sum()}")
print(f"Max start_diameter among ANY multi-record cluster: {multi['start_diameter'].max():,} bp")

print("\n=== Old (misleading) span_bp metric, for comparison ===")
print(f"Multi-record clusters with span_bp_old_metric > 1500bp: "
      f"{(multi['span_bp_old_metric'] > 1500).sum()}")

# Confirm the hypothesis: large span_bp_old_metric clusters are driven by
# individual huge-END records, not by start-position chaining
large_span_small_start_diam = multi[
    (multi["span_bp_old_metric"] > 20000) & (multi["start_diameter"] <= 500)
]
print(f"\nClusters with OLD span>20kb but start_diameter<=500bp "
      f"(= correctly formed, single huge-END record inflating the old metric): "
      f"{len(large_span_small_start_diam)}")

# Check SV type composition of these
large_span_ids = large_span_small_start_diam.index
sv_types_involved = df[df["CLUSTER_ID"].isin(large_span_ids)]["SV_TYPE_NORMALIZED"].value_counts()
print("\nSV types found in these 'falsely flagged' clusters:")
print(sv_types_involved)