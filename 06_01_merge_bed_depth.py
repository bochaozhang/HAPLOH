import pandas as pd
import argparse

# 设置参数解析
parser = argparse.ArgumentParser(description="Merge BED and depth files")
parser.add_argument('--bed_file', required=True, help='Path to the BED file')
parser.add_argument('--depth_file', required=True, help='Path to the depth file')
parser.add_argument('--output_file', required=True, help='Path to the output file')
args = parser.parse_args()

# 读取第一个文件
bed_df = pd.read_csv(args.bed_file, sep='\t', header=None)

# 读取第二个文件
depth_df = pd.read_csv(args.depth_file, sep='\t', header=None)

# 添加列名
bed_df.columns = ['chr', 'region', 'contig', 'pos1', 'pos2', 'index']
depth_df.columns = ['contig', 'index', 'depth']

# 进行匹配
merged_df = pd.merge(bed_df, depth_df, how='left', left_on=['contig', 'index'], right_on=['contig', 'index'])

# 将匹配结果添加到最后一列
merged_df = merged_df[['chr', 'region', 'contig', 'pos1', 'pos2', 'index', 'depth']]

# 保存结果到新的文件
merged_df.to_csv(args.output_file, sep='\t', header=False, index=False)

print(f"匹配结果已保存到 {args.output_file}")



