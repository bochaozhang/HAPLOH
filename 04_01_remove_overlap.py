#!/home/users/bczhang/anaconda3/envs/py3/bin/python

import argparse
import os
import pandas as pd
from Bio import pairwise2

'''
def readCoords(sample,hap):
    """
    Read coords into DataFrame
    """
    folder = "/home/users/xyxu/project/05bczhang/NG/prawnsnew/hifihla_12878/{sample}/HP{hap}_dnadiff".format(sample=sample,hap=hap)
    infile = os.path.join(folder,"out.1coords")
    df = pd.read_csv(infile,sep="\t",names=["S1","E1","S2","E2","LEN 1","LEN 2","% IDY","LEN R","LEN Q","COV R","COV Q","TAG R","TAG Q"])
    return df
'''

def readCoords(workDir, hap):
    folder = os.path.join(workDir, "HP%d_dnadiff" % hap)
    infile = os.path.join(folder, "out.1coords")
    df = pd.read_csv(infile, sep="\t", names=["S1", "E1", "S2", "E2", "LEN 1", "LEN 2", "% IDY", "LEN R", "LEN Q", "COV R", "COV Q", "TAG R", "TAG Q"])
    return df




def filterDirection(df):
    """
    If the majority of the contig hits are in one direction, then filter out those in the other direction
    """
    reverse = [df["S2"][i] > df["E2"][i] for i in range(len(df))] 
    df["REV"] = reverse
    for contig in df["TAG Q"].unique():
        if sum(df["REV"][df["TAG Q"]==contig]) / sum(df["TAG Q"]==contig) > 0.5:
            df = df[~((df["TAG Q"]==contig) & (df["REV"]==False))]
        else:
            df = df[~((df["TAG Q"]==contig) & (df["REV"]==True))]
    return df

def filterCompletSurrounded(df):
    """
    Check if one alignment is completely surrounded in another alignment
    """
    surrounded = [False] * len(df)
    for i in range(len(df)-1):
        for j in range(i+1,len(df)):
            if df["TAG Q"].iloc[i] == df["TAG Q"].iloc[j]:
                if min(df["S2"].iloc[i],df["E2"].iloc[i]) >= min(df["S2"].iloc[j],df["E2"].iloc[j]) and max(df["S2"].iloc[i],df["E2"].iloc[i]) <= max(df["S2"].iloc[j],df["E2"].iloc[j]) and not (min(df["S2"].iloc[i],df["E2"].iloc[i]) == min(df["S2"].iloc[j],df["E2"].iloc[j]) and max(df["S2"].iloc[i],df["E2"].iloc[i]) == max(df["S2"].iloc[j],df["E2"].iloc[j])):
                    surrounded[i] = True
                if min(df["S2"].iloc[i],df["E2"].iloc[i]) <= min(df["S2"].iloc[j],df["E2"].iloc[j]) and max(df["S2"].iloc[i],df["E2"].iloc[i]) >= max(df["S2"].iloc[j],df["E2"].iloc[j]) and not (min(df["S2"].iloc[i],df["E2"].iloc[i]) == min(df["S2"].iloc[j],df["E2"].iloc[j]) and max(df["S2"].iloc[i],df["E2"].iloc[i]) == max(df["S2"].iloc[j],df["E2"].iloc[j])):
                    surrounded[j] = True
    df["surrounded"] = surrounded
    return df
                
def getOverlap(df):
    """
    Check hit pair overlap.
    """
    LEN_THRESHOLD = 1000
    inbetween = lambda x,a,b: True if x in range(min([a,b]),max([a,b])+1) else False
    same_direction = lambda a1,a2,b1,b2: True if (a1-a2) * (b1-b2) >= 0 else False
    overlaps = {"idx1":[], "idx2": [], "pos":[]}
    for contig in df["TAG Q"].unique():
        hits = df[["S1","E1","S2","E2"]][df["TAG Q"]==contig]
        for i in range(len(hits)-1):
            S2_i,E2_i = hits["S2"].iloc[i],hits["E2"].iloc[i]
            for j in range(i+1,len(hits)):
                S2_j,E2_j = hits["S2"].iloc[j],hits["E2"].iloc[j]
                if not same_direction(S2_i,E2_i,S2_j,E2_j):
                    S2_j,E2_j = E2_j,S2_j
                """
                Reserved for reference overlap
                if inbetween(hits["S1"].iloc[j],hits["S1"].iloc[i],hits["E1"].iloc[i]):
                if inbetween(hits["E1"].iloc[j],hits["S1"].iloc[i],hits["E1"].iloc[i]):
                """
                if inbetween(S2_j,S2_i,E2_i):
                    if abs(S2_j - E2_i) > LEN_THRESHOLD:
                        overlaps["idx1"].append(hits.index[i])
                        overlaps["idx2"].append(hits.index[j])
                        overlaps["pos"].append(str(S2_j) + "-" + str(E2_i))
                if inbetween(E2_j,S2_i,E2_i):
                    if abs(S2_i - E2_j) > LEN_THRESHOLD:
                        overlaps["idx1"].append(hits.index[i])
                        overlaps["idx2"].append(hits.index[j])
                        overlaps["pos"].append(str(S2_i) + "-" + str(E2_j))
    overlaps = pd.DataFrame(overlaps)
    return overlaps

def getRefPos(ref_pos,contig_pos,hit_pos):
    """
    Given the hit pair positions on reference and on contig, calculate reference position based on hit postion.
    For example:

    ref    10 11 12 13 14 15 16 17 18 19 20
    contig 3  4  5  6  7  8  9  10 11    12 
    hit 5-7: 10+(5-3)=12, 20-(12-7)=15 --> 12-15

    ref    10 11 12 13 14 15 16 17 18 19 20
    contig 12 11 10 9  8  7  6  5  4     3
    hit 7-5: 10+(12-7)=15, 20-(5-3)=18 --> 15-18
    """
    ref_pos = [int(p) for p in ref_pos]
    contig_pos = [int(p) for p in contig_pos]
    hit_pos = [int(p) for p in hit_pos]
    same_direction = lambda a,b: True if (a[0]-a[1]) * (b[0]-b[1]) >= 0 else False
    if not same_direction(contig_pos,hit_pos):
        hit_pos[0],hit_pos[1] = hit_pos[1],hit_pos[0]
    if contig_pos[0] < contig_pos[1]:
        hit_pos.sort()
        return [ref_pos[0]+(hit_pos[0]-contig_pos[0]), ref_pos[1]-(contig_pos[1]-hit_pos[1])]
    else:
        hit_pos.sort(reverse=True)
        return [ref_pos[0]+(contig_pos[0]-hit_pos[0]), ref_pos[1]-(hit_pos[1]-contig_pos[1])]

def getQueryRegion(ref_pos,contig_pos,hit_pos):
    """
    Reserved for reference overlap
    """
    return 0

def process_file(lines,target_parts,change):
    """
    Credit to Xingyu
    Minor modification by Bochao
    """
    print(target_parts)
    print(change)
    # Adjust the column indices for 0-based indexing
    target_values = [target_parts[11], target_parts[12], target_parts[7], target_parts[8]]
    # 找到以'>'开头且后面的内容为输入行的第12、13、8、9的行
    start_index = None
    for i, line in enumerate(lines):
        if line.startswith('>') and all(value in line for value in target_values):
            start_index = i
            break
    if start_index is None:
        print("Target line not found.")
        return
    # 提取到下一个以'>'开头的行前的内容（若后面没有'>'开头的行就到文件最后一行）
    extracted_lines = []
    for line in lines[start_index + 1:]:
        if line.startswith('>'):
            break
        extracted_lines.append(line.strip())
    # 查找提取行中哪些行有7列
    seven_column_lines = [line for line in extracted_lines if len(line.split()) == 7]
    # 获取第1、2、3、4列为输入行的第1、2、3、4列的行
    first_four_cols_values = target_parts[:4]
    matching_lines_indices = [i for i, line in enumerate(seven_column_lines) if all(line.split()[j] == first_four_cols_values[j] for j in range(4))]
    if not matching_lines_indices:
        print("No matching lines with the first four columns.")
        return
    # 找到第一个匹配的有7列的行和第二个任意有7列的行
    first_match_index = matching_lines_indices[0]
    if len(seven_column_lines) > 1:
        second_match_index = first_match_index + 1 if first_match_index + 1 < len(seven_column_lines) else len(seven_column_lines)
    else:
        second_match_index = len(extracted_lines)
    # 在extracted_lines中的索引
    first_line_index = extracted_lines.index(seven_column_lines[first_match_index])
    if second_match_index < len(seven_column_lines):
        second_line_index = extracted_lines.index(seven_column_lines[second_match_index])
    else:
        second_line_index = len(extracted_lines)
    # 修改行的逻辑
    change_first_value = str(change[0])
    change_second_value = str(change[1])
    # lines[start_index + 1 + first_line_index]中的值
    first_line_parts = lines[start_index + 1 + first_line_index].split()
    # 确定 change_first_value 在哪一列
    if first_line_parts[0] == change_first_value:
        col_index = 0
    elif first_line_parts[1] == change_first_value:
        col_index = 1
    else:
        #############################################
        # added by Bochao
        #############################################
        change[0],change[1] = change[1],change[0]
        change_first_value = str(change[0])
        change_second_value = str(change[1])
        if first_line_parts[0] == change_first_value:
            col_index = 0
        elif first_line_parts[1] == change_first_value:
            col_index = 1
        else:
        #############################################
        # added by Bochao
        #############################################
            print("Change first value not found in the expected columns.")
            return
    if col_index == 0:
        # 修改对应的列值
        first_line_parts[col_index] = change_second_value
        lines[start_index + 1 + first_line_index] = ' '.join(map(str, first_line_parts)) + '\n'
        # 第一列的情况
        dif = abs(change[1] - change[0]) + 1
        sum1 = 0
        sum2 = 0
        line_count = start_index + 2 + first_line_index
        while line_count < start_index + 1 + second_line_index:
            current_value = int(lines[line_count].split()[0])
            if current_value < 0:
                sum1 += abs(current_value) - 1
                sum2 += abs(current_value)
            elif current_value > 0:
                sum1 += current_value
                sum2 += current_value - 1
            elif current_value == 0:
                sum1 += current_value
                sum2 += current_value
            if sum1 >= dif:
                break
            line_count += 1
        third_value = int(first_line_parts[2])
        fourth_value = int(first_line_parts[3])
        if sum1 == dif:
            if third_value < fourth_value:
                third_value = third_value + sum2
            else:
                third_value = third_value - sum2
        elif sum1 > dif:
            modified_line_parts = lines[line_count].split()
            modified_value = int(modified_line_parts[0])
            if third_value < fourth_value:
                if modified_value > 0:
                    third_value = third_value + sum2 - sum1 + dif
                else:
                    third_value = third_value + sum2 - sum1 + dif - 2
            else:
                if modified_value > 0:
                    third_value = third_value - sum2 + sum1 - dif
                else:
                    third_value = third_value - sum2 + sum1 - dif + 2    
        elif sum1 < dif:
            if third_value < fourth_value:
                third_value = third_value + dif - sum1 + sum2 - 1
            else:
                print(sum1)
                print(sum2)
                print(dif)
                print(third_value)
                third_value = third_value - dif + sum1 - sum2 + 1
        if sum1 == dif:
            lines[line_count] = '1\n'
        elif sum1 > dif:
            current_value = int(lines[line_count].split()[0])
            if current_value > 0:
                new_value = sum1 - dif
            else:
                new_value = dif - sum1 - 1
            lines[line_count] = f'{new_value}\n'
        # 修改第三列的值
        first_line_parts[2] = str(third_value)
        lines[start_index + 1 + first_line_index] = ' '.join(map(str, first_line_parts)) + '\n'
        # 更新target_parts
        target_parts[0] = change_second_value
        target_parts[4] = str(int(target_parts[1]) - int(change_second_value) + 1)
        target_parts[2] = str(third_value)
        if int(target_parts[2]) <= int(target_parts[3]):
            target_parts[5] = str(int(target_parts[3]) - int(target_parts[2]) + 1)
        else:
            target_parts[5] = str(int(target_parts[2]) - int(target_parts[3]) + 1)
        # 删除相应的行
        if sum1 >= dif:
            filtered_lines = lines[:start_index + 2 + first_line_index] + lines[line_count:]
        else:
            filtered_lines = lines[:start_index + 2 + first_line_index] + lines[start_index + second_line_index:]
    elif col_index == 1:
        # 第二列的情况
        first_line_parts[col_index] = change_second_value
        lines[start_index + 1 + first_line_index] = ' '.join(map(str, first_line_parts)) + '\n'
        dif = abs(change[1] - int(first_line_parts[0])) + 1
        sum1 = 0
        sum2 = 0
        max_sum1 = 0
        max_sum2 = 0
        line_count = start_index + 2 + first_line_index
        max_line_count = start_index + second_line_index - 1
        max_sum1_line_count = line_count-1
        while line_count <= max_line_count:
            current_value = int(lines[line_count].split()[0])
            if current_value < 0:
                temp_sum1 = sum1 + abs(current_value) - 1
                temp_sum2 = sum2 + abs(current_value)
            else:
                temp_sum1 = sum1 + current_value
                temp_sum2 = sum2 + current_value - 1
            if temp_sum1 <= dif:
                sum1 = temp_sum1
                sum2 = temp_sum2
                max_sum1 = sum1
                max_sum2 = sum2
                max_sum1_line_count = line_count
            if temp_sum1 > dif:
                break
            line_count += 1
        third_value = int(first_line_parts[2])
        fourth_value = int(first_line_parts[3])
        if third_value < fourth_value:
            fourth_value = third_value + dif + max_sum2 - max_sum1 - 1
        else:
            fourth_value = third_value - max_sum2 - dif + max_sum1 + 1
        first_line_parts[3] = str(fourth_value)
        lines[start_index + 1 + first_line_index] = ' '.join(map(str, first_line_parts)) + '\n'
        # 更新target_parts
        target_parts[1] = change_second_value
        target_parts[4] = str(int(change_second_value) - int(target_parts[0]) + 1)
        target_parts[3] = str(fourth_value)
        if int(target_parts[2]) <= int(target_parts[3]):
            target_parts[5] = str(int(target_parts[3]) - int(target_parts[2]) + 1)
        else:
            target_parts[5] = str(int(target_parts[2]) - int(target_parts[3]) + 1)
        # 删除相应的行
        filtered_lines = lines[:max_sum1_line_count + 1] + lines[max_line_count + 1:]
    # 检查并删除连续的以'>'开头的行
    final_lines = []
    previous_was_gt = False
    for line in filtered_lines:
        if line.startswith('>'):
            if previous_was_gt:
                final_lines.pop()  # Remove the previous '>' line
            previous_was_gt = True
        else:
            previous_was_gt = False
        final_lines.append(line)
    return target_parts,final_lines

def delete_line(lines,target_parts):
    """
    Credit to Xingyu
    """
    # Adjust the column indices for 0-based indexing
    target_values = [target_parts[11], target_parts[12], target_parts[7], target_parts[8]]
    # 找到以'>'开头且后面的内容为输入行的第12、13、8、9的行
    start_index = None
    for i, line in enumerate(lines):
        if line.startswith('>') and all(value in line for value in target_values):
            start_index = i
            break
    if start_index is None:
        print("Target line not found.")
        return
    # 提取到下一个以'>'开头的行前的内容（若后面没有'>'开头的行就到文件最后一行）
    extracted_lines = []
    for line in lines[start_index + 1:]:
        if line.startswith('>'):
            break
        extracted_lines.append(line.strip())
    # 查找提取行中哪些行有7列
    seven_column_lines = [line for line in extracted_lines if len(line.split()) == 7]
    # 获取第1、2、3、4列为输入行的第1、2、3、4列的行
    first_four_cols_values = target_parts[:4]
    matching_lines_indices = [i for i, line in enumerate(seven_column_lines) if all(line.split()[j] == first_four_cols_values[j] for j in range(4))]
    if not matching_lines_indices:
        print("No matching lines with the first four columns.")
        return
    # 找到第一个匹配的有7列的行和第二个任意有7列的行
    first_match_index = matching_lines_indices[0]
    if len(seven_column_lines) > 1:
        second_match_index = first_match_index + 1 if first_match_index + 1 < len(seven_column_lines) else len(seven_column_lines)
    else:
        second_match_index = len(extracted_lines)
    # 在extracted_lines中的索引
    first_line_index = extracted_lines.index(seven_column_lines[first_match_index])
    if second_match_index < len(seven_column_lines):
        second_line_index = extracted_lines.index(seven_column_lines[second_match_index])
    else:
        second_line_index = len(extracted_lines)
    # 删除相应的行
    filtered_lines = lines[:start_index + 1] + lines[start_index + 1:start_index + 1 + first_line_index] + lines[start_index + 1 + second_line_index:]
    # 检查并删除连续的以'>'开头的行
    final_lines = []
    previous_was_gt = False
    for line in filtered_lines:
        if line.startswith('>'):
            if previous_was_gt:
                final_lines.pop()  # Remove the previous '>' line
            previous_was_gt = True
        else:
            previous_was_gt = False
        final_lines.append(line)
    return final_lines

def getAlignmentScore(seq1,seq2):
    BIN = 10000
    score = 0
    for i in range(0,min(len(seq1),len(seq2)),BIN):
        sub_seq1 = seq1[i:min(len(seq1),i+BIN)]
        sub_seq2  = seq2[i:min(len(seq2),i+BIN)]
        score += pairwise2.align.globalxx(sub_seq1,sub_seq2,score_only=True)
    return score 

def removeOverlap(lines,df,overlaps,ref_seq,query_seq):
    K = []
    for i in range(len(overlaps)):
        contig = df["TAG Q"].loc[overlaps["idx1"][i]]
        #contig = int(df["TAG Q"].loc[overlaps["idx1"][i]])
        hit_pos = [int(p) for p in overlaps["pos"][i].split("-")]
        if hit_pos[0] < hit_pos[1]:
            seq0 = query_seq[str(contig)][hit_pos[0]-1:hit_pos[1]]
        else:
            seq0 = query_seq[str(contig)][hit_pos[1]-1:hit_pos[0]][::-1]
        ref_pos1 = getRefPos([df["S1"].loc[overlaps["idx1"][i]],df["E1"].loc[overlaps["idx1"][i]]],[df["S2"].loc[overlaps["idx1"][i]],df["E2"].loc[overlaps["idx1"][i]]],hit_pos) 
        ref_pos2 = getRefPos([df["S1"].loc[overlaps["idx2"][i]],df["E1"].loc[overlaps["idx2"][i]]],[df["S2"].loc[overlaps["idx2"][i]],df["E2"].loc[overlaps["idx2"][i]]],hit_pos)
        seq1 = ref_seq[ref_pos1[0]-1:ref_pos1[1]]
        seq2 = ref_seq[ref_pos2[0]-1:ref_pos2[1]]
        print("aligning...")
        score1 = pairwise2.align.globalxx(seq0,seq1,score_only=True)
        score2 = pairwise2.align.globalxx(seq0,seq2,score_only=True)
        remove = 0
        if score1 > score2:
            remove = 2     
        elif score1 < score2:
            remove = 1
        elif score1 == score2:
            if abs(len(seq1)-len(seq0)) > abs(len(seq2)-len(seq0)):
                remove = 1
            elif abs(len(seq1)-len(seq0)) < abs(len(seq2)-len(seq0)):
                remove = 2
            elif abs(len(seq1)-len(seq0)) == abs(len(seq2)-len(seq0)):
                if len(seq1) >= len(seq2):
                    remove = 2
                else:
                    remove = 2
        if remove == 1:
            target_parts = [str(c) for c in df.loc[overlaps["idx1"][i]].tolist()]
            target_parts,lines = process_file(lines,target_parts,ref_pos1) 
            df.loc[overlaps["idx1"][i]] = target_parts
            K.append(ref_pos2)
        elif remove == 2:
            target_parts = [str(c) for c in df.loc[overlaps["idx2"][i]].tolist()]
            target_parts,lines = process_file(lines,target_parts,ref_pos2)
            df.loc[overlaps["idx2"][i]] = target_parts
            K.append(ref_pos1)
        else:
            K.append(0)
            continue
    overlaps["Keep"] = K
    return overlaps,lines



'''
def getRefSeq():
    ref = "/home/users/xyxu/project/05bczhang/NG/prawnsnew/hifihla_12878/10X_hifi/HP1_vcf/hg38.mhc.fa"
    lines = [line.strip() for line in open(ref).readlines()]
    ref_seq = "".join(lines[1:])
    return ref_seq
'''

def getRefSeq(ref):
    lines = [line.strip() for line in open(ref).readlines()]
    ref_seq = "".join(lines[1:])
    return ref_seq



#def getQuerySeq(sample,hap):
#    folder = "/home/users/xyxu/project/05bczhang/NG/prawnsnew/hifihla_12878"
#    fq = os.path.join(folder,sample,"HP%ddedup.50000.filtered.fasta" % hap)
#    lines = [line.strip() for line in open(fq).readlines()]
#    query_seq = {}
#    for i in range(0,len(lines),2):
#        key = int(lines[i].split(" ")[0].replace(">",""))
#        query_seq[key] = lines[i+1]
#    return query_seq




#def getQuerySeq(sample, hap):
#    folder = "/home/users/xyxu/project/05bczhang/NG/prawnsnew/hifihla_12878"
#    fq = os.path.join(folder, sample, "HP%ddedup.50000.filtered.fasta" % hap)
#    lines = [line.strip() for line in open(fq).readlines()]
#    query_seq = {}
#    for i in range(0, len(lines), 2):
#        key = lines[i].split(" ")[0].replace(">", "")
#        query_seq[key] = lines[i+1] if i+1 < len(lines) else ""
#    return query_seq


def getQuerySeq(fq_path):
    #fq = os.path.join(workDir, "HP%ddedup.50000.filtered.fasta" % hap)
    lines = [line.strip() for line in open(fq_path).readlines()]
    query_seq = {}
    for i in range(0, len(lines), 2):
        key = lines[i].split(" ")[0].replace(">", "")
        query_seq[key] = lines[i+1] if i+1 < len(lines) else ""
    return query_seq



def main():
    parser = argparse.ArgumentParser(description="Analyze DNA sequences and coordinates")
    parser.add_argument("workDir", type=str, help="The working directory containing sample data.")
    parser.add_argument("hap1_fq", type=str, help="The path to the hap1 FASTA file.")
    parser.add_argument("hap2_fq", type=str, help="The path to the hap2 FASTA file.")
    #parser.add_argument("sample", type=str, help="The sample ID to analyze.")
    parser.add_argument("ref", type=str, help="The reference sequence file path.")
    args = parser.parse_args()
    workDir = args.workDir
    hap1_fq = args.hap1_fq
    hap2_fq = args.hap2_fq
    query_seq_hap1 = getQuerySeq(hap1_fq)
    query_seq_hap2 = getQuerySeq(hap2_fq)
    #sample = args.sample
    ref_seq = getRefSeq(args.ref)
    #ref_seq = getRefSeq()
    #for hap in [1,2]:
    for hap, query_seq in zip([1, 2], [query_seq_hap1, query_seq_hap2]):
        #query_seq = getQuerySeq(workDir,hap)
        df = readCoords(workDir,hap)
        df = filterCompletSurrounded(df)
        input_file = os.path.join(workDir,"HP%d_dnadiff" % hap,"out.1delta")
        lines = open(input_file).readlines()
        #df = filterDirection(df)
        for i in range(len(df)):
            if df["surrounded"].iloc[i]:
                target_parts = [str(c) for c in df.loc[i].tolist()] 
                lines = delete_line(lines,target_parts)
        df = df[~df["surrounded"]]
        overlaps = getOverlap(df)  
        overlaps,lines = removeOverlap(lines,df,overlaps,ref_seq,query_seq)
        print("writing to output...") 
        output_overlap = os.path.join(workDir,"HP%d_dnadiff" % hap,"overlaps.txt")
        overlaps.to_csv(output_overlap, index=False)
        output_file = os.path.join(workDir,"HP%d_dnadiff" % hap,"out_noOverlap.1delta")
        with open(output_file,'w') as o:
            o.writelines(line for line in lines)


if __name__ == "__main__":
    main()
