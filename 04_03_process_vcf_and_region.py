import pandas as pd
import sys
from pathlib import Path

# ===============================
# 参数
# ===============================
sample = sys.argv[1]
workDir = Path(sys.argv[2])

# ===============================
# Step 1: together.vcf → together_all.vcf
# ===============================
together_vcf = workDir / "together.vcf"
together_all_vcf = workDir / "together_all.vcf"

pon_single_df = pd.read_csv(together_vcf, sep="\t", header=None)

# 添加 VCF 固定列
pon_single_df.insert(len(pon_single_df.columns), "QUAL", ".")
pon_single_df.insert(len(pon_single_df.columns), "FILTER", "PASS")
pon_single_df.insert(len(pon_single_df.columns), "INFO", "AF=0.5")
pon_single_df.insert(len(pon_single_df.columns), "FORMAT", "OP")

# 调整列顺序
pon_single_df = pon_single_df.iloc[:, list(range(0, 5)) + list(range(6, 10)) + [5]]

# VCF header
header = (
    "##fileformat=VCFv4.1\n"
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    "##INFO=<ID=AF,Number=1,Type=Float,Description=\"assembled haploid frequency\">\n"
    "##INFO=<ID=OP,Number=.,Type=String,Description=\"Original Position\">\n"
    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
)

with open(together_all_vcf, "w") as f:
    f.write(header)
    pon_single_df.to_csv(f, sep="\t", index=False, header=False)

print(f"[OK] Generated {together_all_vcf}")

# ===============================
# Step 2: region 过滤 VCF
# ===============================
region_file = workDir / f"{sample}.filtered.region.hg38.txt"
together_all_filter_vcf = workDir / "together_all.filter.vcf"

regions = pd.read_csv(region_file, sep="\t", header=None)

with open(together_all_vcf, "rt") as f:
    lines = f.readlines()

header_lines = [l for l in lines if l.startswith("#")]
data_lines = [l for l in lines if not l.startswith("#")]

vcf_data = pd.DataFrame([l.strip().split("\t") for l in data_lines])
vcf_data[1] = vcf_data[1].astype(int)

mask = pd.Series(False, index=vcf_data.index)

for _, row in regions.iterrows():
    start, end = row[0], row[1]
    mask |= (vcf_data[1] >= start) & (vcf_data[1] <= end)

filtered_vcf_data = vcf_data[mask]

with open(together_all_filter_vcf, "wt") as f:
    for line in header_lines:
        f.write(line)
    filtered_vcf_data.to_csv(f, sep="\t", header=False, index=False)

print(f"[OK] Generated {together_all_filter_vcf}")

# ===============================
# Step 3: 生成 hc_region.txt
# ===============================
hc_region_file = workDir / "hc_region.txt"

df = pd.read_csv(
    region_file,
    sep="\t",
    header=None,
    usecols=[0, 1, 2, 3, 12]
)

df["chr"] = "chr6"
df = df[["chr", 0, 1, 12, 2, 3]]

df.to_csv(hc_region_file, sep="\t", header=False, index=False)

print(f"[OK] Generated {hc_region_file}")
print("[DONE] All steps finished successfully.")


