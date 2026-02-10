import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
from functools import reduce
# Function to get purity
def getPurity(sample,purity_csv):
    df = pd.read_csv(purity_csv)
    purity = float(df["tumor_purity"][df["sample"] == sample])
    return purity

# Function to plot regions from a file
def plot_region_lines(file_path, color, region_label):
    if os.path.getsize(file_path) > 0:  # Check if the file is not empty
        region_data = pd.read_csv(file_path, sep='\t', header=None)
        for _, row in region_data.iterrows():
            plt.axvline(x=row[0], color=color, linestyle='--')
            plt.axvline(x=row[1], color=color, linestyle='--')
        plt.plot([], [], color=color, linestyle='--', label=region_label)

# Main function to run the analysis
def main(sample, window_size,smooth_size,purity_csv,regions_to_analyze,pos_range):
    # Get purity
    purity = getPurity(sample,purity_csv)
    # Load CSV file
    file_path = f'{sample}.csv'
    df = pd.read_csv(file_path, sep='\t')
    df.dropna(inplace=True)
    # Define regions to analyze
    #regions_to_analyze = [(28903952, 30915858),(30915859,31177951), (31177952, 33268517)]
    # Initialize M values dictionary
    M_values = {}
    # Calculate M for each region
    for start, end in regions_to_analyze:
        region_df = df[(df['pos'] >= start) & (df['pos'] <= end)]
        normal_sum = region_df['normal_cov'].sum()
        tumor_sum = region_df['tumor_cov'].sum()
        M_values[(start, end)] = normal_sum / tumor_sum if tumor_sum != 0 else 0
    # Filter the DataFrame to include only the relevant regions
    #filtered_df = df[(df['pos'] >= regions_to_analyze[0][0]) & (df['pos'] <= regions_to_analyze[0][1]) |
    #                 (df['pos'] >= regions_to_analyze[1][0]) & (df['pos'] <= regions_to_analyze[1][1]) |
    #                 (df['pos'] >= regions_to_analyze[2][0]) & (df['pos'] <= regions_to_analyze[2][1])]
    filtered_df = df[reduce(lambda x, y: x | y, [(df['pos'] >= start) & (df['pos'] <= end) for start, end in regions_to_analyze])]
    filtered_df = filtered_df[filtered_df['normal_cov'] != 0]
    # Calculate ploidy for each position
    ploidy_values = []
    for index, row in filtered_df.iterrows():
        pos = row['pos']
        tumor_cov = row['tumor_cov']
        normal_cov = row['normal_cov']
        # Determine which region the current position falls into
        for (start, end), M in M_values.items():
            if start <= pos <= end:
                # Calculate ploidy
                ploidy = (2 * tumor_cov / normal_cov * M - 2 + 2 * purity) / purity
                ploidy_values.append((pos, ploidy))
                break
    # Convert ploidy values to DataFrame
    ploidy_df = pd.DataFrame(ploidy_values, columns=['pos', 'ploidy'])
    # Load smoothed ploidy file
    file_path = f'{sample}_Ploidy_sum_{int(window_size/500)}kbp.txt'
    ploidy_df_smooth = pd.read_csv(file_path, sep='\t')
    #pos_range = [28953952, 33268517]
    ploidy_filtered = ploidy_df[(ploidy_df['pos'] >= pos_range[0]) & (ploidy_df['pos'] <= pos_range[1])]
    ploidy_smooth_filtered = ploidy_df_smooth[(ploidy_df_smooth['pos'] >= pos_range[0]) & (ploidy_df_smooth['pos'] <= pos_range[1])]
    # Apply window averaging
    ploidy_filtered = ploidy_filtered.copy()
    ploidy_filtered['window'] = (ploidy_filtered['pos'] // smooth_size) * smooth_size
    ploidy_smooth_filtered = ploidy_smooth_filtered.copy()
    ploidy_smooth_filtered['window'] = (ploidy_smooth_filtered['pos'] // smooth_size) * smooth_size
    # Calculate averages
    ploidy_avg = ploidy_filtered.groupby('window').mean().reset_index()
    ploidy_smooth_avg = ploidy_smooth_filtered.groupby('window').mean().reset_index()
    # Merge the averages into one DataFrame
    merged_avg_df = pd.merge(ploidy_avg, ploidy_smooth_avg, on='window', suffixes=('_x', '_y'))
    # Plot ploidy values
    plt.figure(figsize=(15, 6))
    plt.scatter(merged_avg_df['pos_x'], merged_avg_df['ploidy_x'], label='Ploidy X (avg)', color='lightblue', s=1)
    plt.plot(ploidy_smooth_filtered['pos'], ploidy_smooth_filtered['ploidy'], label='Ploidy Y (smoothed avg)', linestyle='--', color='gray')
    # Plot deletion and amplification regions
    low_file_path = f'{sample}.low.sum.{int(window_size/500)}kbp.txt'
    high_file_path = f'{sample}.high.sum.{int(window_size/500)}kbp.txt'
    plot_region_lines(low_file_path, 'purple', 'Deletion Region')
    plot_region_lines(high_file_path, 'green', 'Amplification Region')
    # Set plot limits and labels
    plt.ylim(-2, 4)
    plt.xlabel('Position (pos)')
    plt.ylabel('Ploidy Value')
    plt.title('Ploidy Values with Averaging')
    # Display legend outside plot
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left')
    # Save and show plot
    plt.savefig(f'{sample}_plo_avg_loess.png')
    plt.show()
    # Process coverage data for normal and tumor
    cov_df = df[(df['pos'] >= pos_range[0]) & (df['pos'] <= pos_range[1])]
    cov_df = cov_df.copy()
    cov_df['window'] = (cov_df['pos'] // smooth_size) * smooth_size
    cov_df = cov_df[['window', 'pos', 'normal_cov', 'tumor_cov', 'normal_cov_smooth', 'tumor_cov_smooth']].groupby('window').mean().reset_index()
    # Calculate dynamic ylim based on max and 90th percentile of 'normal_cov'
    normal_cov_90th = cov_df['normal_cov'].quantile(0.9)
    tumor_cov_90th = cov_df['tumor_cov'].quantile(0.9)
    cov_ylim_upper = int(max(1.5 * tumor_cov_90th, 1.5 * normal_cov_90th))
    # Plot coverage values
    plt.figure(figsize=(15, 6))
    plt.scatter(cov_df['pos'], cov_df['normal_cov'], label='Normal Cov (avg)', color='cornflowerblue', s=1)
    plt.plot(df['pos'], df['normal_cov_smooth'], label='Normal Cov (smoothed)', linestyle='--', color='darkblue')
    plt.scatter(cov_df['pos'], cov_df['tumor_cov'], label='Tumor Cov (avg)', color='lightcoral', s=1)
    plt.plot(df['pos'], df['tumor_cov_smooth'], label='Tumor Cov (smoothed)', linestyle='--', color='darkred')
    # Plot deletion and amplification regions in the coverage plot as well
    plot_region_lines(low_file_path, 'purple', 'Deletion Region')
    plot_region_lines(high_file_path, 'green', 'Amplification Region')
    # Set dynamic ylim
    plt.ylim(-50, cov_ylim_upper)
    # Set coverage plot limits and labels
    plt.xlabel('Position (pos)')
    plt.ylabel('Cov Value')
    plt.title('Cov Values with Averaging')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    # Save and show plot
    plt.savefig(f'{sample}_cov_avg_sum_{int(window_size/500)}kbp.png')
    plt.show()

# Setup argparse for command line input
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze ploidy and coverage data.')
    parser.add_argument('--sample', type=str, required=True, help='Sample identifier')
    parser.add_argument('--window_size', type=int, default=5000, help='Window size for smoothing (default: 5000)')
    parser.add_argument('--smooth_size', type=int, default=2000, help='Window size for averaging')
    parser.add_argument('--purity_csv', type=str, required=True, help='Path to the tumor purity CSV file')
    parser.add_argument('--start', type=int, default=28903952, help='Start position for range (default: 28903952)')
    parser.add_argument('--end', type=int, default=33268518, help='End position for range (default: 33268518)')
    #parser.add_argument('--regions_to_analyze', type=str, default="[(28903952,33268517)]",
    #                    help='Regions to analyze as a list of tuples (default: [(28903952,33268517)])')
    parser.add_argument('--regions_to_analyze', type=str,
                        default="[(%(start)s, %(end)s)]",
                        help='Regions to analyze as a list of tuples, default is [(start, end)]')
    args = parser.parse_args()
    regions_to_analyze = eval(args.regions_to_analyze)
    pos_range = [args.start,args.end]
    main(args.sample, args.window_size,args.smooth_size,args.purity_csv,regions_to_analyze,pos_range)
