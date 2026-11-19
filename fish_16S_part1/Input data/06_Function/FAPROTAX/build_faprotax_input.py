import csv
import re
import os
import sys

base = "/public/home/liruojing/fish_16S_merge/analysis_A_global/eukaryota_filtered"
outdir = os.path.join(base, "FAPROTAX")

table_file = os.path.join(outdir, "feature_table.tsv")
tax_file = os.path.join(base, "taxonomy_no_eukaryota.tsv")
output_file = os.path.join(outdir, "faprotax_input.tsv")

# 读取 taxonomy
taxonomy = {}
with open(tax_file, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        feature_id = row["Feature ID"].strip()
        taxon = row["Taxon"].strip()
        taxonomy[feature_id] = taxon

# 读取 ASV 丰度表
with open(table_file, "r", encoding="utf-8-sig") as f:
    next(f)  # 跳过 BIOM 注释行
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    samples = header[1:]
    rows = list(reader)

output_rows = []
matched = 0
unmatched = []
domain_count = {"Bacteria": 0, "Archaea": 0, "Other": 0}

for row in rows:
    if not row:
        continue

    feature_id = row[0].strip()
    counts = row[1:]

    if feature_id not in taxonomy:
        unmatched.append(feature_id)
        continue

    taxon = taxonomy[feature_id]
    levels = [x.strip() for x in taxon.split(";") if x.strip()]

    if not levels:
        continue

    domain = re.sub(r"^[a-z]__", "", levels[0]).strip()

    if domain not in ("Bacteria", "Archaea"):
        domain_count["Other"] += 1
        continue

    domain_count[domain] += 1

    clean_taxonomy = [
        re.sub(r"^[a-z]__", "", x).strip()
        for x in levels
    ]
    clean_taxonomy = [x for x in clean_taxonomy if x]

    if not clean_taxonomy:
        continue

    output_rows.append(["|".join(clean_taxonomy)] + counts)
    matched += 1

if unmatched:
    print("错误：存在未匹配的 ASV，未生成最终输入文件。")
    print("未匹配数量：", len(unmatched))
    print("示例：", unmatched[:10])
    sys.exit(1)

with open(output_file, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["Taxonomy"] + samples)
    writer.writerows(output_rows)

print("FAPROTAX 输入文件生成成功")
print("丰度表 ASV 数量：", len(rows))
print("成功匹配并保留的细菌/古菌 ASV：", matched)
print("Bacteria：", domain_count["Bacteria"])
print("Archaea：", domain_count["Archaea"])
print("其他分类：", domain_count["Other"])
print("样本数量：", len(samples))
print("未匹配 ASV：", len(unmatched))
print("输出文件：", output_file)
