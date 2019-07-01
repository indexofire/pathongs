# MLST

---

传统的 MLST 分析方法采用对样本进行 PCR 扩增后进行 Sanger 测序、拼接、比对后，获得等位基因序列，根据碱基差异分配序号，再根据相应的序号组合确定 ST 型別。常用的工具有bigsDB，Bionumerics等。

随着 NGS 技术应用与微生物基因组测序，越来越多的细菌基因组被测序生成草图。几乎所有的数据都可以覆盖 MLST 涉及的管家基因，因此你可以直接通过基因组上管家基因的序列来分析获得其ST型別。

[MLST](https://github.com/tseemann/mlst) 是由 [@Torsten Seemann](https://twitter.com/torstenseemann) 用 perl 开发的实现从fasta序列获取MLST型別的工具。

## 1. 安装

```bash
# mlst 工具在 bioconda 有安装包，可以直接通过 conda 安装
$ conda create -n mlst
$ conda activate mlst
(mlst)$ conda install mlst
```

## 2. 使用

MLST 支持 `.fa`,`.fasta`,`.gb`,`.gbk`,`.fna`以及这些格式的gz压缩包格式的文件。不加参数可以对所有物种数据库进行扫描，加上参数`--scheme`可以对指定物种进行扫描。

```bash
# 直接扫描 fasta 数据
(mlst)$ mlst mygenome.fasta

# 支持 multiple fasta 格式
(mlst)$ mlst scaffolds.fasta

# 支持多文件格式
(mlst)$ mlst data/*

# 查看可以扫描的物种
(mlst)$ mlst --list
# 使用 listera monocytogenes 物种数据库扫描
(mlst)$ mlst --scheme lmonocytogenes data/*
```

## 3. 应用举例

对 NCBI genome 数据库中 bacillus cereus 的 Assembly 进行扫描获取 MLST 数据。

```bash
# edirect 抓取 FTP 链接，awk补全下载路径，用 wget 工具批量下载（设置相应参数避免进入下载假死状态）
$ esearch -db assembly -query "Bacillus cereus[ORGN] AND latest[SB]" | \
> efetch -format docsum | \
> xtract -pattern DocumentSummary -element FtpPath_RefSeq | \
> awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}' | \
> xargs wget --no-passive-ftp --limit-rate=200k

# mlst 扫描
(mlst)$ mlst --scheme bcereus *.fna.gz > ../bcereus-mlst-result.txt
(mlst)$ cat ../bcereus-mlst-result.txt
```
