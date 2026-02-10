import sys
import pandas as pd

def extract_data(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                split_line = line.strip().split('\t')
                data.append([split_line[0], int(split_line[1]), split_line[3], split_line[4], split_line[9]])
    return data

def read_bed(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'region', 'col3', 'pos', 'col5', 'col6'])

def get_pos(num, data):
    try:
        row = data.loc[num]
        return f"{row['col3']}-{row['col6']}"
    except KeyError:
        return float('nan')

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: python script.py hp1_data hp2_data hp1_bed hp2_bed output_file")
        sys.exit(1)

    hp1_data_file, hp2_data_file, hp1_bed_file, hp2_bed_file, output_file = sys.argv[1:6]

    hp1_data = extract_data(hp1_data_file)
    hp2_data = extract_data(hp2_data_file)

    hp1_df = pd.DataFrame(hp1_data, columns=['chrom', 'hg38_pos', 'hg38.nt', 'hp1.nt', 'hp1_pos'])
    hp2_df = pd.DataFrame(hp2_data, columns=['chrom', 'hg38_pos', 'hg38.nt', 'hp2.nt', 'hp2_pos'])

    merged_df = pd.merge(hp1_df, hp2_df, on=['chrom', 'hg38_pos', 'hg38.nt'], how='outer')

    hp1_bed = read_bed(hp1_bed_file)
    hp2_bed = read_bed(hp2_bed_file)

    hp1_bed.set_index('pos', inplace=True)
    hp2_bed.set_index('pos', inplace=True)

    for index, row in merged_df.iterrows():
        if pd.isna(row['hp1.nt']) and pd.isna(row['hp1_pos']):
            pos = get_pos(row['hg38_pos'], hp1_bed)
            if not pd.isna(pos):
                merged_df.at[index, 'hp1.nt'] = row['hg38.nt']
                merged_df.at[index, 'hp1_pos'] = pos
        if pd.isna(row['hp2.nt']) and pd.isna(row['hp2_pos']):
            pos = get_pos(row['hg38_pos'], hp2_bed)
            if not pd.isna(pos):
                merged_df.at[index, 'hp2.nt'] = row['hg38.nt']
                merged_df.at[index, 'hp2_pos'] = pos

    to_drop = []
    for index, row in merged_df.iterrows():
        #if row['hp1.nt'] == row['hp2.nt']:
        #    to_drop.append(index)
        if pd.isna(row['hp1.nt']) or pd.isna(row['hp2.nt']):
            to_drop.append(index)

    merged_df.drop(to_drop, inplace=True)
    merged_df = merged_df[['chrom', 'hg38_pos', 'hp1_pos', 'hp2_pos', 'hg38.nt', 'hp1.nt', 'hp2.nt']]
    merged_df.to_csv(output_file, sep='\t', index=False)

