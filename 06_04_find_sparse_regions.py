import pandas as pd
import argparse

def find_non_overlapping_regions(vcf_file, start, end, threshold, min_region_size, output_file):
    # 读取VCF文件
    vcf_df = pd.read_csv(vcf_file, sep='\t')

    # 获取变异位点的位置信息，并确保它们在给定范围内
    positions = vcf_df['hg38_pos'].values
    positions = positions[(positions >= start) & (positions <= end)]

    # 初始化结果列表
    regions = []
    current_start = start

    while current_start < end:
        # 初始region的结束位置
        current_end = current_start + min_region_size

        # 扩展当前region，直到包含的变异位点数目达到threshold
        in_region = positions[(positions >= current_start) & (positions < current_end)]

        while len(in_region) <= threshold and current_end < end:
            current_end += 1
            in_region = positions[(positions >= current_start) & (positions < current_end)]

        # 添加满足条件的region到结果列表
        if current_end - current_start > min_region_size:
            regions.append((current_start, current_end - 1))

        # 移动到下一个region
        current_start = current_end

    # 合并相邻的区域（如果它们的间隔为1）
    merged_regions = []
    for region in regions:
        if merged_regions and merged_regions[-1][1] + 1 == region[0]:
            # 如果当前region与前一个region相邻，则合并它们
            merged_regions[-1] = (merged_regions[-1][0], region[1])
        else:
            merged_regions.append(region)

    # 保存结果到文件
    with open(output_file, 'w') as f_out:
        for region_start, region_end in merged_regions:
            f_out.write(f"{region_start}\t{region_end}\n")

    print(f"结果已保存到 {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find non-overlapping regions with variant counts below a threshold and merge adjacent regions")
    parser.add_argument('--vcf_file', required=True, help='Path to the VCF file')
    parser.add_argument('--start', type=int, required=True, help='Start position of the search range')
    parser.add_argument('--end', type=int, required=True, help='End position of the search range')
    parser.add_argument('--threshold', type=int, default=2, help='Maximum number of variants allowed in a region')
    parser.add_argument('--min_region_size', type=int, default=1000, help='Minimum size of each region')
    parser.add_argument('--output_file', required=True, help='Path to the output file')
    args = parser.parse_args()

    find_non_overlapping_regions(args.vcf_file, args.start, args.end, args.threshold, args.min_region_size, args.output_file)

