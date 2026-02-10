import pandas as pd
import matplotlib.pyplot as plt
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description="Analyze ploidy data and extract significant high/low regions.")
parser.add_argument('--sample', type=str, required=True, help='Sample ID to be analyzed, e.g., 1645442')
parser.add_argument('--start', type=int, default=28903952, help='Start position for range (default: 28903952)')
parser.add_argument('--end', type=int, default=33268518, help='End position for range (default: 33268518)')
#parser.add_argument('--regions_to_analyze', type=str, default="[(28903952,33268517)]",
#                    help='Regions to analyze as a list of tuples (default: [(28903952,33268517)])')
parser.add_argument('--regions_to_analyze', type=str,
                    default="[(%(start)s, %(end)s)]",
                    help='Regions to analyze as a list of tuples, default is [(start, end)]')


parser.add_argument('--min_length', type=int, default=25000, help='min deletion length (default: 25000)')
parser.add_argument('--min_valid_points', type=int, default=12500, help='min valid points (default: 12500)')
parser.add_argument('--window_size', type=int, default=5000, help='Window size for moving average')


args = parser.parse_args()

# 获取样本ID
sample = args.sample

# 读取文件
file_path = f'{sample}_Ploidy_sum_{int(args.window_size/500)}kbp.txt'

# 假设文件是tab分隔的，读取数据
df = pd.read_csv(file_path, sep='\t')

# 分析的区域
#regions_to_analyze = [(28903952, 30915858), (31177952, 33268517)]
#regions_to_analyze = [(28903952, 30915858),(30915859,31177951), (31177952, 33268517)]
regions_to_analyze = eval(args.regions_to_analyze)
# 初始化列表来存储显著的区域
regions_of_interest = []

# 设置最小标准
min_length = args.min_length  # 最小长度
min_valid_points = args.min_valid_points  # 最少有效点数

# 分析每个指定的区域
for start, end in regions_to_analyze:
    # 获取整个指定区域的数据
    segment = df[(df['pos'] >= start) & (df['pos'] <= end)].copy()
    if segment.empty:
        continue
    # 计算整个区域的平均覆盖率和标准差
    mean_cov = segment['ploidy'].mean()
    std_cov = segment['ploidy'].std()
    # 定义显著偏差的阈值：上下两个方向（超过均值的两倍标准差）
    #upper_threshold = mean_cov + 2 * std_cov
    #lower_threshold = mean_cov - 2 * std_cov
    lower_threshold = min(1, mean_cov - 2 * std_cov)
    upper_threshold = max(3, mean_cov + 2 * std_cov)
    # 初始化变量来存储当前的子区域
    current_region = []
    extending = False  # 标记是否在扩展区域
    # 遍历每个数据点
    for i, row in segment.iterrows():
        if current_region:
            # 判断是否要停止扩展
            if (row['ploidy'] <= upper_threshold) and (row['ploidy'] >= lower_threshold):
                # 如果当前点不显著，停止扩展，记录显著区域
                if len(current_region) >= min_valid_points and (current_region[-1]['pos'] - current_region[0]['pos']) >= min_length:
                    sub_segment = pd.DataFrame(current_region)
                    significant_regions_high = sub_segment[sub_segment['ploidy'] > upper_threshold]
                    significant_regions_low = sub_segment[sub_segment['ploidy'] < lower_threshold]
                    if not significant_regions_high.empty:
                        regions_of_interest.append(('high', sub_segment['pos'].iloc[0], sub_segment['pos'].iloc[-1], significant_regions_high))
                    if not significant_regions_low.empty:
                        regions_of_interest.append(('low', sub_segment['pos'].iloc[0], sub_segment['pos'].iloc[-1], significant_regions_low))
                # 重置 current_region 开始新的积累
                current_region = []
                extending = False
        # 如果是显著点或者正在扩展中，继续扩展
        if (row['ploidy'] > upper_threshold) or (row['ploidy'] < lower_threshold) or extending:
            current_region.append(row)
            extending = True  # 继续扩展
    # 处理最后一个子区域，如果还有剩余
    if current_region and len(current_region) >= min_valid_points:
        sub_segment = pd.DataFrame(current_region)
        # 判断是否有显著高于或低于的区域
        significant_regions_high = sub_segment[sub_segment['ploidy'] > upper_threshold]
        significant_regions_low = sub_segment[sub_segment['ploidy'] < lower_threshold]
        if not significant_regions_high.empty:
            regions_of_interest.append(('high', sub_segment['pos'].iloc[0], sub_segment['pos'].iloc[-1], significant_regions_high))
        if not significant_regions_low.empty:
            regions_of_interest.append(('low', sub_segment['pos'].iloc[0], sub_segment['pos'].iloc[-1], significant_regions_low))

# 将high和low区域分别保存到文件
with open(f"{sample}.high.sum.{int(args.window_size/500)}kbp.txt", "w") as high_file, open(f"{sample}.low.sum.{int(args.window_size/500)}kbp.txt", "w") as low_file:
    for region in regions_of_interest:
        level = region[0]   # 'high' or 'low'
        start_pos = region[1]
        end_pos = region[2]
        if level == 'high':
            high_file.write(f"{start_pos}\t{end_pos}\n")
        elif level == 'low':
            low_file.write(f"{start_pos}\t{end_pos}\n")

# 生成图像
plt.figure(figsize=(12, 6))

# Plot all positions in gray
plt.scatter(df['pos'], df['ploidy'], color='gray', s=10, label='All Positions', alpha=0.5)

# Variables to control the legend
high_plotted = False
low_plotted = False

# Overlay significant regions with different colors
for region in regions_of_interest:
    level = region[0]   # 'high' or 'low'
    data = region[3]    # The DataFrame containing the data
    # Debugging: Check the content of each region
    print(f"Region level: {level}, Number of points: {len(data)}")
    color = 'red' if level == 'high' else 'blue'  # Red for high, blue for low
    # Check if data is not empty and plot
    if not data.empty:
        if level == 'high':
            plt.scatter(data['pos'], data['ploidy'], color=color, s=10, 
                        label='Significant High' if not high_plotted else "", alpha=0.8)
            high_plotted = True
        elif level == 'low':
            plt.scatter(data['pos'], data['ploidy'], color=color, s=10, 
                        label='Significant Low' if not low_plotted else "", alpha=0.8)
            low_plotted = True

plt.xlabel('Position')
plt.ylabel('Ploidy')
plt.title(f'Ploidy with Significant High and Low Regions for Sample {sample}')
plt.legend()
plt.grid(True)
plt.savefig(f'{sample}_Significant_Ploidy_Regions_sum_{int(args.window_size/500)}kbp.png')
plt.show()

