import argparse
import pandas as pd

# =========================
# 工具函数
# =========================
def load_hp_data(path, sample):
    hp1_file = f"{path}/HP1_dnadiff_nooverlap/{sample}.HP1.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed"
    hp2_file = f"{path}/HP2_dnadiff_nooverlap/{sample}.HP2.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed"
    df_hp1 = pd.read_csv(hp1_file, sep="\t", header=None)
    df_hp2 = pd.read_csv(hp2_file, sep="\t", header=None)
    return df_hp1, df_hp2


def get_pos_both(df_hp1, df_hp2, pos):
    r1 = df_hp1[df_hp1[3] == pos]
    r2 = df_hp2[df_hp2[3] == pos]
    if not r1.empty and not r2.empty:
        return (
            r1.iloc[0, 2], r1.iloc[0, 5],
            r2.iloc[0, 2], r2.iloc[0, 5],
        )
    return None


def get_pos_single(df, pos):
    row = df[df[3] == pos]
    if not row.empty:
        return row.iloc[0, 2], row.iloc[0, 5]
    return None, None


# =========================
# Step 1：寻找 overlap region
# =========================
def find_overlap_regions(path, sample):
    file1 = f"{path}/HP1_dnadiff_nooverlap/hc_region.txt"
    file2 = f"{path}/HP2_dnadiff_nooverlap/hc_region.txt"
    df1 = pd.read_csv(file1, sep="\t", header=None)
    df2 = pd.read_csv(file2, sep="\t", header=None)
    # 保证 start < end
    df1.loc[df1[1] >= df1[2], [1, 2]] = df1.loc[df1[1] >= df1[2], [2, 1]].values
    df2.loc[df2[1] >= df2[2], [1, 2]] = df2.loc[df2[1] >= df2[2], [2, 1]].values
    df1 = df1.sort_values(by=1)
    df2 = df2.sort_values(by=1)
    df_hp1, df_hp2 = load_hp_data(path, sample)
    overlap_regions = []
    for _, row in df1.iterrows():
        overlaps = df2[
            (df2[0] == row[0]) &
            (df2[1] <= row[2]) &
            (df2[2] >= row[1])
        ]
        for _, o in overlaps.iterrows():
            start = int(max(row[1], o[1]))
            end = int(min(row[2], o[2]))
            overlap_regions.append((row[0], start, end))
    # 边界收缩（保证 start/end 在 HP1 & HP2 都有 pos）
    final_regions = []
    for chrom, start, end in overlap_regions:
        start0 = start
        for p in range(start, end + 1):
            if get_pos_both(df_hp1, df_hp2, p):
                start = p
                break
        for p in range(end, start0 - 1, -1):
            if get_pos_both(df_hp1, df_hp2, p):
                end = p
                break
        if start <= end:
            final_regions.append((chrom, start, end))
    return pd.DataFrame(final_regions, columns=[0, 1, 2])


# =========================
# Step 2：补充 pos 信息
# =========================
def add_pos_info(df_regions, df_hp1, df_hp2):
    df = df_regions.copy()
    df[[3, 4]] = df.apply(lambda r: get_pos_single(df_hp1, r[1]), axis=1, result_type="expand")
    df[[5, 6]] = df.apply(lambda r: get_pos_single(df_hp1, r[2]), axis=1, result_type="expand")
    df[[7, 8]] = df.apply(lambda r: get_pos_single(df_hp2, r[1]), axis=1, result_type="expand")
    df[[9,10]] = df.apply(lambda r: get_pos_single(df_hp2, r[2]), axis=1, result_type="expand")
    df[[3,5,7,9]] = df[[3,5,7,9]].astype(str)
    df[[4,6,8,10]] = df[[4,6,8,10]].astype(int, errors="ignore")
    return df


# =========================
# main
# =========================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Base path")
    parser.add_argument("sample", help="Sample name")
    args = parser.parse_args()
    # Step 1
    overlap_df = find_overlap_regions(args.path, args.sample)
    overlap_file = f"{args.path}/{args.sample}.nooverlap.overlap.region.txt"
    overlap_df.to_csv(overlap_file, sep="\t", index=False, header=False)
    # Step 2
    df_hp1, df_hp2 = load_hp_data(args.path, args.sample)
    final_df = add_pos_info(overlap_df, df_hp1, df_hp2)
    out_file = f"{args.path}/{args.sample}.nooverlap.overlap.region.pos.txt"
    final_df.to_csv(out_file, sep="\t", index=False, header=False)
    print(f"[OK] Output written:")
    print(f"  - {overlap_file}")
    print(f"  - {out_file}")


if __name__ == "__main__":
    main()

