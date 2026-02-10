import pandas as pd
import re
import sys
from pathlib import Path

def get_pos(value, df):
    result = df[df[4] == value][5]
    if not result.empty:
        return result.values[0]
    else:
        return None

def float_to_int(x):
    if isinstance(x, float) and x.is_integer():
        return int(x)
    return x

def process_sample(sample, workDir):
    workDir = Path(workDir)
    coords_file = workDir / f'{sample}.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed'
    df_coords = pd.read_csv(coords_file, sep='\t', header=None)
    out_file = workDir / 'out_noOverlap.1coords'
    df_out = pd.read_csv(out_file, sep='\t', header=None)
    for i in range(1, len(df_out)):
        if df_out.iloc[i-1, 1] >= df_out.iloc[i, 0]:
            temp = df_out.iloc[i, 0]
            original_second_column_value = df_out.iloc[i-1, 1]
            while True:
                temp -= 1
                pos = get_pos(temp, df_coords)
                if pos is not None:
                    df_out.iloc[i-1, 1] = temp
                    df_out.iloc[i-1, 3] = pos
                    break
            temp = original_second_column_value
            while True:
                temp += 1
                pos = get_pos(temp, df_coords)
                if pos is not None:
                    df_out.iloc[i, 0] = temp
                    df_out.iloc[i, 2] = pos
                    break
    region_file = workDir / f'{sample}.region.txt'
    df_out.to_csv(region_file, sep='\t', header=False, index=False)
    df = df_out.copy()
    df_dropped_list = []
    indices_to_drop = []
    for i in range(len(df)):
        if i > 0 and i < len(df) - 1:
            if df.iloc[i, -1] != df.iloc[i-1, -1] and df.iloc[i, -1] != df.iloc[i+1, -1]:
                df_dropped_list.append(df.iloc[i])
                indices_to_drop.append(df.index[i])
        elif i == 0:
            if df.iloc[i, -1] != df.iloc[i+1, -1]:
                df_dropped_list.append(df.iloc[i])
                indices_to_drop.append(df.index[i])
        elif i == len(df) - 1:
            if df.iloc[i, -1] != df.iloc[i-1, -1]:
                df_dropped_list.append(df.iloc[i])
                indices_to_drop.append(df.index[i])
    if df_dropped_list:
        df_dropped = pd.concat(df_dropped_list, axis=1).T
    else:
        df_dropped = pd.DataFrame()
    if not df_dropped.empty:
        df_dropped = df_dropped[df_dropped.iloc[:, 4] < 1500]
        df_dropped = df_dropped.applymap(float_to_int)
    df_dropped.to_csv(workDir / 'potential_error_region.txt', sep='\t', header=False, index=False)
    df_dropped_set = set([tuple(row) for row in df_dropped.values])
    df_filtered = df[~df.apply(tuple, axis=1).isin(df_dropped_set)]
    df_filtered.to_csv(workDir / f'{sample}.filtered.region.txt', sep='\t', header=False, index=False)
    print(f"[DONE] Sample={sample} processed. Region and filtered.region files saved.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_regions.py <sample> <workDir>")
        sys.exit(1)
    sample = sys.argv[1]
    workDir = sys.argv[2]
    process_sample(sample, workDir)




