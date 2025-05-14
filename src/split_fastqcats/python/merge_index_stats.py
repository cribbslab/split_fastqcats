import pandas as pd
import os
import argparse

def merge_stats_by_index(input_dir, output_file):
    # Get the sample name by stripping the '.index_stats.csv' from the output filename
    sample_name = os.path.basename(output_file).replace('.index_stats.csv', '')

    # Initialize a list to store the stats files
    stats_files = []

    # Walk through the input directory and its subdirectories
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            # Only include files that start with the sample name and end with '.stats.csv'
            if file.startswith(sample_name) and file.endswith(".stats.csv"):
                stats_files.append(os.path.join(root, file))

    # Initialize a dictionary to store merged index counts
    merged_index_counts = {}
    merged_metrics = {}

    # Process each stats file
    for stats_file in stats_files:
        # Read the CSV file
        df = pd.read_csv(stats_file, header=None)

        # First table (First 7 rows)
        # Assuming the first 7 rows contain metric data
        metric_df = df.iloc[:7]
        metric_df.columns = metric_df.iloc[0]  # Set the first row as column headers
        metric_df = metric_df[1:]
        # Second table (Next 13 rows)
        # Assuming the next 13 rows contain index counts
        index_df = df.iloc[7:20]
        index_df.columns = index_df.iloc[0]  # Set the first row as column headers
        index_df = index_df[1:]

        # Process metrics

        # Aggregate metrics (sum values)
        for _, row in metric_df.iterrows():
            metric = row['Metric']
            count = int(row['Value'])  # Convert to float
            merged_metrics[metric] = merged_metrics.get(metric, 0) + count


        # Process index counts (the second table)
        for _, row in index_df.iterrows():
            index = row['Index']
            count = int(row['SegmentCount'])  # Ensure this is an integer
            merged_index_counts[index] = merged_index_counts.get(index, 0) + count

    # Convert the merged metrics to DataFrame
    metrics_df = pd.DataFrame(list(merged_metrics.items()), columns=["Metric", "Value"])
    # Convert the merged index counts to DataFrame
    index_df = pd.DataFrame(list(merged_index_counts.items()), columns=["Index", "SegmentCount"])
    
    # Write both metrics and index counts to the output file (tab-separated)
    with open(output_file, "w") as f:
        # Write the metrics first
        metrics_df.to_csv(f, sep=",", index=False)
        f.write("\n")  # Add a newline between the two tables
        # Write the index counts
        index_df.to_csv(f, sep=",", index=False)
    # Write the merged statistics to the output file (tab-separated)
    #index_df.to_csv(output_file, sep=",", index=False)

def merge_stats_by_index(input_dir, output_file):
    import os
    import pandas as pd
    from collections import Counter

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
