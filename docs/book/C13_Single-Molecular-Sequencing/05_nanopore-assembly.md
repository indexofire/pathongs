# Nanopore Reads 组装

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## minimap2 + miniasm

- minimap2
- miniasm
- Minion_QC

### 建立环境

```bash
$ conda create -n minasm
$ conda activate minasm
(minasm)$ conda install minimap2 miniasm Minion_QC
```

### 组装流程

```bash
# QC
(minasm)$ Rscript MinIONQC.R -i data/reads.fastq -o qc/result

# 过滤段片段，切除5'和3'端低质量序列
(minasm)$ NanoFilt -l 1000 --headcrop 30 -q 9 < reads.fastq > reads_trimmed.fastq

# 比对序列
(minasm)$ minimap -x ava-ont reads_trimmed.fastq | gzip -1 > minimap.paf.gz

# 组装序列
(minasm)$ miniasm -f reads_trimmed.fastq minimap.paf.gz > miniasm.gfa
(minasm)$ awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fasta
(minasm)$ assembly-stats miniasm.fasta

# mummer 比较2个基因组共线性
(mummer)$ dnadiff -p result ref.fasta miniasm.fasta
# 访问 www.assemblytics.com，用可视化工具查看比较接俄国
```

## Flye 组装

```bash
# 组装
(flye)$ flye --nano-raw reads_trimmed.fastq --genome-size 5m --out-dir flye_result

# mummer 比较2个基因组
(mummer)$ dnadiff -p result ref.fasta flye_result/assembly.fasta
```

## Reference

1. https://timkahlke.github.io/LongRead_tutorials
