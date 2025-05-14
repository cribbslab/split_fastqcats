import os
import pandas as pd
import argparse
from collections import Counter


def merge_stats(input_dir, output_file):
    """
    Merge all individual stats CSV files into a single stats file by summing the values for each metric.
    The final output will be saved in the output_file path.
    """
    metrics_dict = {
        'total_reads': 0,
        'processed_reads': 0,
        'total_segments': 0,
        'full_length_segments': 0,
        'lowqual_segments': 0,
        'binned_reads_0': 0,
	'binned_reads_1':0
    }
    
    segment_hit_counts = Counter()

    sample_name = os.path.basename(output_file).replace('.stats.csv', '')
    
    
    

    # Use os.walk to traverse directories recursively
    for root, dirs, files in os.walk(input_dir):
        for stats_file in files:
            if stats_file.startswith(sample_name) and stats_file.endswith(".stats.csv"):
                stats_filepath = os.path.join(root, stats_file)
                
                # Read the stats file into a DataFrame
                df = pd.read_csv(stats_filepath)
                
                # Sum up the metrics from the file and update the total in the metrics_dict
                for metric in metrics_dict:
                    if metric in df['Metric'].values:
                        metric_value = df[df['Metric'] == metric]['Value'].values[0]
                        metrics_dict[metric] += metric_value
    
    # Convert the merged results into a DataFrame
    merged_stats = pd.DataFrame(list(metrics_dict.items()), columns=["Metric", "Summed Value"])
    
    # Save the merged stats to a CSV file
    merged_stats.to_csv(output_file, index=False)


def merge_stats(input_dir, output_file):
    metrics_dict = {
        'total_reads': 0,
        'processed_reads': 0,
        'total_segments': 0,
        'full_length_segments': 0,
        'lowqual_segments': 0,
        'binned_reads_0': 0,
        'binned_reads_1': 0
    }

    segment_hit_counts = Counter()
    primer_hit_counts = Counter()
    sample_name = os.path.basename(output_file).replace('.stats.csv', '')

    for root, dirs, files in os.walk(input_dir):
        for stats_file in files:
            if stats_file.startswith(sample_name) and stats_file.endswith(".stats.csv"):
                stats_filepath = os.path.join(root, stats_file)
                
                with open(stats_filepath, 'r') as f:
                    lines = f.readlines()

                section = "metrics"

                for line in lines:
                    line = line.strip()
                    if not line:
                        if section == "metrics":
                            section = "segments"
                        elif section == "segments":
                            section = "primers"
                        continue

                    if line.startswith("Metric") or line.startswith("Segment Hits") or line.startswith("Primer Hits"):
                        continue

                    try:
                        key, value = line.split(",")
                        value = int(value)

                        if section == "metrics" and key in metrics_dict:
                            metrics_dict[key] += value
                        elif section == "segments":
                            segment_hit_counts[int(key)] += value
                        elif section == "primers":
                            primer_hit_counts[int(key)] += value
                    except ValueError:
                        continue  # skip malformed lines

    # Write merged output
    with open(output_file, 'w', newline='') as f:
        f.write("Metric,Summed Value\n")
        for metric, value in metrics_dict.items():
            f.write(f"{metric},{value}\n")

        f.write("\nSegment Hits,Read Count\n")
        for hit_count in sorted(segment_hit_counts):
            f.write(f"{hit_count},{segment_hit_counts[hit_count]}\n")

        f.write("\nPrimer Hits,Read Count\n")
        for hit_count in sorted(primer_hit_counts):
            f.write(f"{hit_count},{primer_hit_counts[hit_count]}\n")


def main():
    parser = argparse.ArgumentParser(description="Merge individual stats CSV files into a single merged file by summing the values.")
    parser.add_argument("--input-dir", required=True, help="Directory containing individual stats CSV files (including subdirectories).")
    parser.add_argument("--output", required=True, help="Path to the output merged stats file.")
     
    args = parser.parse_args()

    # Run the merging function
    merge_stats(args.input_dir, args.output)

if __name__ == "__main__":
    main()
