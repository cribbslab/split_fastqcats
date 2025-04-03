import os
import argparse
import pandas as pd

def merge_stats_by_index(input_dir, output_file):
    """
    Merge all individual stats CSV files from input_dir into a single CSV file.
    The final output will be saved in output_file.
    """
    stats_files = []

    # Find all stats files
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".stats.csv"):
                stats_files.append(os.path.join(root, file))

    # Storage for combined metrics
    merged_metrics = {}
    merged_index_counts = {}

    for stats_file in stats_files:
        df = pd.read_csv(stats_file, sep="\t")  # Assuming tab-separated values

        # Separate metrics and index counts
        for _, row in df.iterrows():
            key = row.iloc[0]  # First column contains the metric/index
            value = row.iloc[1]  # Second column contains the value

            if key.startswith("total_") or key in ["processed_reads", "lowqual_segments", "binned_reads"]:
                merged_metrics[key] = merged_metrics.get(key, 0) + value
            else:  # These are index counts
                merged_index_counts[key] = merged_index_counts.get(key, 0) + value

    # Convert to DataFrame
    metrics_df = pd.DataFrame(list(merged_metrics.items()), columns=["Metric", "Value"])
    index_df = pd.DataFrame(list(merged_index_counts.items()), columns=["Index", "SegmentCount"])

    # Write to output file
    output_path = os.path.join(output_file, "merged_index_stats.csv")
    with open(output_path, "w") as f:
        metrics_df.to_csv(f, sep="\t", index=False)
        f.write("\n")  # Separate sections
        index_df.to_csv(f, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Merge index stats from all chunks into one file.")
    parser.add_argument("--input-dir", required=True, help="Directory containing individual index stats CSV files.")
    parser.add_argument("--output-file", required=True, help="Directory to save the merged stats file.")

    args = parser.parse_args()

    merge_stats_by_index(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()
