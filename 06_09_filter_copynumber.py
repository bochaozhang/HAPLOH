import argparse
import pandas as pd

# 定义命令行参数
def parse_args():
    parser = argparse.ArgumentParser(description="Filter copy number data based on region exclusion.")
    parser.add_argument("--id", required=True, help="Sample ID (e.g., 1640487)")
    parser.add_argument("--region_file", required=True, help="Path to the region file (e.g., low.sum.10kbp.txt)")
    parser.add_argument("--copynumber_file", required=True, help="Path to the copy number file (e.g., copynumber_loess_150bp_sum.csv)")
    parser.add_argument("--output_file", required=True, help="Path to the output file (filtered result)")
    return parser.parse_args()

# 主函数
def main():
    args = parse_args()
    sample_id = args.id
    region_file = args.region_file
    copynumber_file = args.copynumber_file
    output_file = args.output_file
    # 读取区域文件，第1、2列构成一个 region
    region_df = pd.read_csv(region_file, sep="\t", header=None, names=["start", "end"])
    # 读取拷贝数文件
    copynumber_df = pd.read_csv(copynumber_file)
    # 确保 hg38_pos 列为数字类型
    copynumber_df["hg38_pos"] = pd.to_numeric(copynumber_df["hg38_pos"], errors="coerce")
    # 去除 hg38_pos 在 region 范围内的行
    filtered_df = copynumber_df[~copynumber_df["hg38_pos"].apply(
        lambda pos: any((pos >= row.start and pos <= row.end) for row in region_df.itertuples())
    )]
    # 输出到文件
    filtered_df.to_csv(output_file, index=False)
    print(f"Filtered results saved to {output_file}")

# 运行主函数
if __name__ == "__main__":
    main()


