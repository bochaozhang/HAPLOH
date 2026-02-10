import pandas as pd
import argparse

# 设置参数解析
parser = argparse.ArgumentParser(description="Merge two BED files by selected columns and sum a specific column")
parser.add_argument('--file_hap1', required=True, help='Path to the first BED file (hap1)')
parser.add_argument('--file_hap2', required=True, help='Path to the second BED file (hap2)')
parser.add_argument('--output_file', required=True, help='Path to the output file')
args = parser.parse_args()

# 读取两个文件，只选择需要的列
hap1_df = pd.read_csv(args.file_hap1, sep='\t', header=None, usecols=[0, 3, 6])
hap2_df = pd.read_csv(args.file_hap2, sep='\t', header=None, usecols=[0, 3, 6])

# 为每个DataFrame命名列
hap1_df.columns = ['chrom', 'pos', 'hap1_value']
hap2_df.columns = ['chrom', 'pos', 'hap2_value']

# 合并两个DataFrame，只保留chrom和pos都存在的行
merged_df = pd.merge(hap1_df, hap2_df, on=['chrom', 'pos'])

# 将hap1_value和hap2_value相加
merged_df['sum_value'] = merged_df['hap1_value'] + merged_df['hap2_value']

# 选择需要的列并保存结果
output_df = merged_df[['chrom', 'pos', 'sum_value']]
output_df.to_csv(args.output_file, sep='\t', header=False, index=False)

print(f"The merged result has been saved to {args.output_file}")

