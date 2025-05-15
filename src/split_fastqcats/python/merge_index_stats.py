import pandas as pd
import os
import argparse
from collections import Counter



def merge_stats_by_index(input_dir, output_file):

    sample_name = os.path.basename(output_file).replace('.index_stats.csv', '')

    stats_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.startswith(sample_name) and file.endswith(".stats.csv"):
                stats_files.append(os.path.join(root, file))

    merged_metrics = {}
    merged_index_counts = Counter()
    merged_index_hits = Counter()

    for stats_file in stats_files:
        with open(stats_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]  # Strip and remove blank lines

        # Identify section headers
        section_indices = {
            'metrics': None,
            'index_counts': None,
            'index_hits': None
        }

        for i, line in enumerate(lines):
            if line.startswith("Metric,Value"):
                section_indices['metrics'] = i
            elif line.startswith("Index,SegmentCount"):
                section_indices['index_counts'] = i
            elif line.startswith("Index hits,Read count"):
                section_indices['index_hits'] = i

        # Metrics section
        metrics_start = section_indices['metrics'] + 1
        metrics_end = section_indices['index_counts']
        for line in lines[metrics_start:metrics_end]:
            metric, value = line.split(",")
            merged_metrics[metric] = merged_metrics.get(metric, 0) + int(value)

        # Index counts section
        index_start = section_indices['index_counts'] + 1
        index_end = section_indices['index_hits']
        for line in lines[index_start:index_end]:
            idx, count = line.split(",")
            merged_index_counts[idx] += int(count)

        # Index hits section
        hits_start = section_indices['index_hits'] + 1
        for line in lines[hits_start:]:
            hit, count = line.split(",")
            merged_index_hits[int(hit)] += int(count)

    # Output
    with open(output_file, "w") as f:
        # Metrics
        f.write("Metric,Value\n")
        for k, v in merged_metrics.items():
            f.write(f"{k},{v}\n")
        f.write("\n")

        # Index counts
        f.write("Index,SegmentCount\n")
        for k, v in merged_index_counts.items():
            f.write(f"{k},{v}\n")
        f.write("\n")

        # Index hits
        f.write("Index hits,Read Count\n")
        for k in sorted(merged_index_hits):
            f.write(f"{k},{merged_index_hits[k]}\n")



def main():
    parser = argparse.ArgumentParser(description="Merge index stats from all chunks into one file.")
    parser.add_argument("--input-dir", required=True, help="Directory containing individual index stats CSV files.")
    parser.add_argument("--output-file", required=True, help="Output file to save the merged stats.")

    args = parser.parse_args()
    merge_stats_by_index(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()
