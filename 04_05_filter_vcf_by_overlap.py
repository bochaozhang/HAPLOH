import pandas as pd
import sys
from pathlib import Path

# ===============================
# 参数
# ===============================
sample = sys.argv[1]
workdir = Path(sys.argv[2])

# ===============================
# 输入文件路径
# ===============================
region_file = workdir.parent / f"{sample}.nooverlap.overlap.region.pos.txt"
vcf_file = workdir / "together_all.vcf"

# ===============================
# 输出文件路径
# ===============================
out_vcf = workdir / "together_all.overlap.filter.vcf"

# ===============================
# 读取 region 文件
# ===============================
regions = pd.read_csv(region_file, sep="\t", header=None)

# ===============================
# 读取 VCF（保留 header）
# ===============================
with open(vcf_file, "rt") as f:
    lines = f.readlines()

header_lines = [line for line in lines if line.startswith("#")]
data_lines = [line for line in lines if not line.startswith("#")]

vcf_data = pd.DataFrame([line.strip().split("\t") for line in data_lines])
vcf_data[1] = vcf_data[1].astype(int)

# ===============================
# 构建 mask
# ===============================
mask = pd.Series(False, index=vcf_data.index)

for _, row in regions.iterrows():
    start, end = row[1], row[2]
    mask |= (vcf_data[1] >= start) & (vcf_data[1] <= end)

filtered_vcf_data = vcf_data[mask]

# ===============================
# 写出结果
# ===============================
with open(out_vcf, "wt") as f:
    for line in header_lines:
        f.write(line)
    filtered_vcf_data.to_csv(f, sep="\t", header=False, index=False)

print(f"[OK] Written: {out_vcf}")

