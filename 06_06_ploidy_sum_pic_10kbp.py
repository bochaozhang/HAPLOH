import pandas as pd
import matplotlib.pyplot as plt
import argparse
import ast
# Function to get purity based on sample ID
def getPurity(sample,purity_csv):
    #df = pd.read_csv("/home/users/bczhang/HLALOH/tumor_purity.csv")
    df = pd.read_csv(purity_csv)
    purity = float(df["tumor_purity"][df["sample"] == sample])
    return purity


def calculate_smooth_coverage_within_regions(df, regions_to_analyze, window_size=5000):
    """
    计算指定区域内的肿瘤覆盖率和正常覆盖率的平滑值。
    """
    smoothed_values = []
    for start, end in regions_to_analyze:
        # 筛选出当前区域的数据
        region_data = df[(df['pos'] >= start) & (df['pos'] <= end)].copy()
        positions = region_data['pos'].values
        tumor_cov = region_data['tumor_cov'].values
        normal_cov = region_data['normal_cov'].values
        
        # 初始化平滑值列表
        tumor_cov_smooth = []
        normal_cov_smooth = []
        
        # 逐位置计算平滑值
        for pos in positions:
            mask = (positions >= pos - window_size) & (positions <= pos + window_size)
            tumor_cov_smooth.append(tumor_cov[mask].sum())
            normal_cov_smooth.append(normal_cov[mask].sum())
        
        # 将结果存储到当前区域数据中
        region_data['tumor_cov_smooth'] = tumor_cov_smooth
        region_data['normal_cov_smooth'] = normal_cov_smooth
        smoothed_values.append(region_data)
    
    # 合并所有区域数据
    return pd.concat(smoothed_values)

def analyze_sample(sample_id,purity_csv,regions_to_analyze,window_size=5000):
    # Get purity
    purity = getPurity(sample_id,purity_csv)
    # Load CSV file dynamically based on sample_id
    file_path = f'{sample_id}.csv'
    df = pd.read_csv(file_path, sep='\t')
    df.dropna(inplace=True)
    
    # Define regions to analyze
    #regions_to_analyze = [(28903952, 30915858), (30915859, 31177951), (31177952, 33268517)]
    
    # Precompute M values for each region
    M_values = {}
    for start, end in regions_to_analyze:
        region_df = df[(df['pos'] >= start) & (df['pos'] <= end)]
        normal_sum = region_df['normal_cov'].sum()
        tumor_sum = region_df['tumor_cov'].sum()
        M_values[(start, end)] = normal_sum / tumor_sum if tumor_sum != 0 else 0
    
    # Smooth tumor and normal coverage within regions
    df = calculate_smooth_coverage_within_regions(df, regions_to_analyze, window_size=5000)
    
    # Drop positions outside regions to ensure only specified regions are analyzed
    df = df[df['tumor_cov_smooth'].notna()]
    
    # Assign region index
    for i, (start, end) in enumerate(regions_to_analyze):
        df.loc[(df['pos'] >= start) & (df['pos'] <= end), 'region'] = i

    # Calculate ploidy
    df['ploidy'] = df.apply(
        lambda row: (2 * row['tumor_cov_smooth'] / row['normal_cov_smooth'] * 
                     M_values[regions_to_analyze[int(row['region'])]] - 2 + 2 * purity) / purity
        if row['normal_cov_smooth'] > 0 else np.nan, axis=1
    )
    
    # Save results
    ploidy_df = df[['pos', 'ploidy']].dropna()
    ploidy_df.to_csv(f'{sample_id}_Ploidy_sum_{int(window_size/500)}kbp.txt', sep='\t', index=False)
    
    # Plot pos vs. ploidy
    plt.figure(figsize=(10, 6))
    plt.plot(ploidy_df['pos'], ploidy_df['ploidy'], marker='o', linestyle='-', color='b')
    plt.title(f'Ploidy vs. Position for Sample {sample_id}')
    plt.xlabel('Position')
    plt.ylabel('Ploidy')
    plt.grid(True)
    plt.savefig(f'{sample_id}_Ploidy_sum_{int(window_size/500)}kbp.png')
    plt.show()


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Analyze ploidy for a given sample.")
    parser.add_argument("--sample_id", type=str, required=True, help="The sample ID to analyze.")
    parser.add_argument('--purity_csv', type=str, required=True, help='Path to the tumor purity CSV file')
    parser.add_argument('--window_size', type=int, default=5000, help='Window size for smoothing (default: 5000)')
    # Parse arguments
    parser.add_argument('--start', type=int, default=28903952, help='Start position for range (default: 28903952)')
    parser.add_argument('--end', type=int, default=33268518, help='End position for range (default: 33268518)')
    #parser.add_argument('--regions_to_analyze', type=str, default="[(28903952,33268517)]",
#                    help='Regions to analyze as a list of tuples (default: [(28903952,33268517)])')
    parser.add_argument('--regions_to_analyze', type=str,
                        default="[(%(start)s, %(end)s)]",
                        help='Regions to analyze as a list of tuples, default is [(start, end)]')
    args = parser.parse_args()
    regions_to_analyze = eval(args.regions_to_analyze)
    #print(args.regions_to_analyze)
    #print(type(regions_to_analyze))  # 期望是 <class 'list'>
    # Analyze the sample provided by the user
    analyze_sample(args.sample_id,args.purity_csv,regions_to_analyze,args.window_size)


