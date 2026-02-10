import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np

# 设置参数解析
parser = argparse.ArgumentParser(description="Merge BED files, calculate coverage, and plot smoothed coverage")
parser.add_argument('--start', type=int, default=28903952, help='Start position for range (default: 28903952)')
parser.add_argument('--end', type=int, default=33268518, help='End position for range (default: 33268518)')
#parser.add_argument('--regions_to_analyze', type=str, default="[(28903952,33268517)]", 
#                    help='Regions to analyze as a list of tuples (default: [(28903952,33268517)])')
parser.add_argument('--regions_to_analyze', type=str, 
                    default="[(%(start)s, %(end)s)]", 
                    help='Regions to analyze as a list of tuples, default is [(start, end)]')
parser.add_argument('--file_hp1', required=True, help='Path to the HP1 BED file')
parser.add_argument('--file_hp2', required=True, help='Path to the HP2 BED file')
parser.add_argument('--output_file', required=True, help='Path to the output file')
parser.add_argument('--cov_1_png', required=True, help='Path to save the first coverage plot')
parser.add_argument('--cov_miss_png', required=True, help='Path to save the missing coverage plot')
parser.add_argument('--cov_log_png', required=True, help='Path to save the log ratio plot')
parser.add_argument('--hp1_col', required=True, help='Column name for hap1 normal coverage')
parser.add_argument('--hp2_col', required=True, help='Column name for hap1 tumor coverage')
parser.add_argument('--average_size', type=int, default=25000, help='Window size for moving average')
args = parser.parse_args()

# 读取文件
hp1_df = pd.read_csv(args.file_hp1, sep='\t', header=None, names=['chrom', 'pos', 'hap1_value' ], usecols=[0, 1, 2])
hp2_df = pd.read_csv(args.file_hp2, sep='\t', header=None, names=['chrom', 'pos', 'hap2_value'], usecols=[0, 1, 2])

# 重命名列
hp1_df.columns = ['chrom', 'pos', args.hp1_col]
hp2_df.columns = ['chrom', 'pos', args.hp2_col]

# 合并数据框
merged_df = pd.merge(hp1_df, hp2_df, on=['chrom', 'pos'], how='outer').dropna(subset=['chrom', 'pos'])
merged_df[args.hp2_col] = merged_df[args.hp2_col].astype(float)

# 生成从28903952到33268517的完整pos序列
#all_positions = pd.DataFrame({'pos': range(28903952, 33268518)})
all_positions = pd.DataFrame({'pos': range(args.start, args.end)})
# 合并生成的pos序列和已有数据
merged_df = pd.merge(all_positions, merged_df, on='pos', how='left')

# 移动平均平滑函数
def moving_average(df, average_size):
    rolling_window = df.rolling(window=average_size, min_periods=1, center=True)
    return rolling_window.mean()

# 计算移动平均，平滑覆盖度

def smooth_within_regions(df, pos_col, value_col, regions, average_size):
    df[f'{value_col}_smooth'] = pd.NA  # Initialize smooth column with NaNs
    for region in regions:
        # Select the rows that belong to the current region
        region_mask = (df[pos_col] >= region[0]) & (df[pos_col] <= region[1])
        region_df = df.loc[region_mask, value_col]
        # Apply moving average only to the rows within the current region
        df.loc[region_mask, f'{value_col}_smooth'] = moving_average(region_df, average_size * 2)
    return df

# 假设 'pos' 列代表位置，'hp1_col' 和 'hp2_col' 是需要平滑的列

#regions_to_analyze = [(28903952, 30915858), (30915859, 31177951), (31177952, 33268517)]
regions_to_analyze = eval(args.regions_to_analyze)
# 对 hp1_col 和 hp2_col 分别进行区间内的平滑
merged_df = smooth_within_regions(merged_df, 'pos', args.hp1_col, regions_to_analyze, args.average_size)
merged_df = smooth_within_regions(merged_df, 'pos', args.hp2_col, regions_to_analyze, args.average_size)

#merged_df[f'{args.hp1_col}_smooth'] = moving_average(merged_df[args.hp1_col], args.window_size * 2)
#merged_df[f'{args.hp2_col}_smooth'] = moving_average(merged_df[args.hp2_col], args.window_size * 2)
# 将原本没有的 pos 的平滑值设为 NA
merged_df.loc[merged_df[args.hp1_col].isna(), f'{args.hp1_col}_smooth'] = np.nan
merged_df.loc[merged_df[args.hp2_col].isna(), f'{args.hp2_col}_smooth'] = np.nan
merged_df[f'{args.hp1_col}_smooth'].replace(0, np.nan, inplace=True)
merged_df[f'{args.hp2_col}_smooth'].replace(0, np.nan, inplace=True)
# 在计算 log_ratio 之前，针对只存在有效数值的行
# 确保平滑列是 float 类型，并且没有无效数据
merged_df[f'{args.hp1_col}_smooth'] = pd.to_numeric(merged_df[f'{args.hp1_col}_smooth'], errors='coerce')
merged_df[f'{args.hp2_col}_smooth'] = pd.to_numeric(merged_df[f'{args.hp2_col}_smooth'], errors='coerce')

# 创建有效掩码，检查两列是否有有效的数值（非 NaN 和 非 inf）
valid_mask = merged_df[f'{args.hp1_col}_smooth'].notna() & merged_df[f'{args.hp2_col}_smooth'].notna()
valid_mask &= np.isfinite(merged_df[f'{args.hp1_col}_smooth']) & np.isfinite(merged_df[f'{args.hp2_col}_smooth'])

# 检查是否有除零情况，并过滤掉这些行
#valid_mask &= (merged_df[f'{args.hp1_col}_smooth'] != 0) & (merged_df[f'{args.hp2_col}_smooth'] != 0)

merged_df.loc[valid_mask, 'log_ratio'] = np.where(
    merged_df.loc[valid_mask,f'{args.hp1_col}_smooth'] == 0,
    np.nan,
    np.log2(merged_df.loc[valid_mask,f'{args.hp2_col}_smooth'] / merged_df.loc[valid_mask,f'{args.hp1_col}_smooth'])
)
# 对有效数据进行 log_ratio 计算
merged_df.loc[valid_mask, 'log_ratio'] = np.log2(merged_df.loc[valid_mask, f'{args.hp2_col}_smooth'] / merged_df.loc[valid_mask, f'{args.hp1_col}_smooth'])

# 对 log_ratio 进行平滑处理
merged_df['log_ratio_smooth'] = merged_df['log_ratio'].rolling(window=args.average_size * 2, min_periods=1, center=True).mean()

#merged_df['log_ratio'] = np.log2(merged_df[f'{args.hp2_col}_smooth'] / merged_df[f'{args.hp1_col}_smooth'])
#merged_df['log_ratio_smooth'] = merged_df['log_ratio'].rolling(window=args.window_size * 2, min_periods=1, center=True).mean()
#merged_df = merged_df.dropna()
merged_df = merged_df.dropna(subset=['chrom'])


# 绘制并保存第一个覆盖图
plt.figure(figsize=(12, 6))
plt.plot(merged_df['pos'], merged_df[f'{args.hp1_col}_smooth'], label=f'{args.hp1_col}_smooth', linestyle='-')
plt.plot(merged_df['pos'], merged_df[f'{args.hp2_col}_smooth'], label=f'{args.hp2_col}_smooth', linestyle='-')
plt.xlabel('Position')
plt.ylabel('Smoothed Coverage')
plt.title(f'Smoothed Coverage vs. Position (Using {args.average_size * 2} bp Window)')
plt.legend()
plt.grid(True)
plt.savefig(args.cov_1_png)

# 对缺失值进行重新索引并绘制图表
merged_df1 = merged_df.set_index('pos').reindex(range(merged_df['pos'].min(), merged_df['pos'].max() + 1)).reset_index()
plt.figure(figsize=(12, 6))
plt.plot(merged_df1['pos'], merged_df1[f'{args.hp1_col}_smooth'], label=f'{args.hp1_col}_smooth', linestyle='-')
plt.plot(merged_df1['pos'], merged_df1[f'{args.hp2_col}_smooth'], label=f'{args.hp2_col}_smooth', linestyle='-')
plt.xlabel('Position')
plt.ylabel('Smoothed Coverage')
plt.title(f'Smoothed Coverage vs. Position (Using {args.average_size * 2} bp Window)')
plt.legend()
plt.grid(True)
plt.savefig(args.cov_miss_png)

# 计算比值并取对数
#merged_df1['log_ratio'] = np.log2(merged_df1[f'{args.hp2_col}_smooth'] / merged_df1[f'{args.hp1_col}_smooth'])

# 对 log_ratio 进行再次平滑
#merged_df1['log_ratio_smooth'] = merged_df1['log_ratio'].rolling(window=args.window_size * 2, min_periods=1, center=True).mean()
plt.figure(figsize=(12, 6))
plt.plot(merged_df1['pos'], merged_df1['log_ratio_smooth'], label='log_ratio_smooth', linestyle='-')
plt.xlabel('Position')
plt.ylabel('Log2 Ratio Smoothed')
plt.title(f'Smoothed Log2 Ratio ({args.hp2_col} / {args.hp1_col}) vs. Position (Using {args.average_size * 2} bp Window)')
plt.legend()
plt.grid(True)
min_val = merged_df1['log_ratio_smooth'].min()
max_val = merged_df1['log_ratio_smooth'].max()

if min_val >= -2 and max_val <= 2:
    plt.ylim(-2, 2)
plt.savefig(args.cov_log_png)

# 保存匹配结果到输出文件
merged_df.to_csv(args.output_file, sep='\t', index=False)

print(f"The smoothed plot has been saved as {args.cov_1_png}, {args.cov_miss_png}, and {args.cov_log_png}")


