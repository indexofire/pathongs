# ARIBA

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    [Ariba](https://github.com/sanger-pathogens/ariba) 工具是由 Sanger 测序中心病原微生物开发组开发的用于基因组数据扫描的工具。主要用于耐药基因的扫描，也可以用来对特定基因扫描（如MLST管家基因）。

## 1. 安装

版本：

- python: 3.6.12
- ariba: 2.14.6
- bowtie2: 2.3.5
- cd-hit: 4.8.1
- mummer: 3.23

```bash
# 新建虚拟环境
$ conda create -n ariba
$ conda activate ariba

# 安装ariba及其依赖包
(ariba)$ conda install ariba
```

## 2. 使用

**获得数据库**

```bash
# 定义数据库路径
$ export DBS=your-database-path

# 下载耐药和毒力基因数据库文件
(ariba)$ mkdir $DBS/ariba
(ariba)$ for i in argannot card megares plasmidfinder resfinder srst2_argannot vfdb_core vfdb_full virulencefinder
> do ariba getref $i $DBS/ariba/$i
> done
# 格式化数据库文件生成ariba可用的数据库
# 这里以card为例，其他数据库进行相应格式化
(ariba)$ ariba prepareref -f $DBS/ariba/card.fa -m $DBS/ariba/card.tsv card

# 下载PubMLST数据库
(ariba)$ mkdir $DBS/ariba/pubmlst
(ariba)$ ariba pubmlstspecies > $DBS/ariba/pubmlst/species.txt
# 批量下载所有物种pubmlst数据信息
(ariba)$ IFS=$'\n'; for i in `cat $DBS/ariba/pubmlst/species.txt`; do ariba pubmlstget $i $DBS/ariba/pubmlst/$i; done
# PubMLST数据在用pubmlstget过程中已经格式化完成，可以直接使用。
```

**分析数据**

```bash
# 分析 reads mlst
(ariba)$ ariba run "$DBS/ariba/pubmlst/Salmonella enterica/ref_db" S1_1.fastq.gz S1_2.fastq.gz S1_output

# 分析耐药数据
(ariba)$ ariba run "$DBS/ariba/card" S1_1.fastq.gz S1_2.fastq.gz S1_output
```

**生成多个样本报告**

```bash
(ariba)$ for i in `ls *_1.fastq.gz`; do ariba run $DBS/ariba/card $i ${i%_1.fastq.gz}_2.fastq.gz ${i%_1.fastq.gz}_output; done
(ariba)$ ariba summary result *_output/report*.tsv
```

**构建自定义数据库**

```bash
# 准备数据库序列文件和meta数据

# 格式化序列和meta数据生成数据库
(ariba)$ ariba prepareref -f myseq.fa -m meta.tsv mydata
(ariba)$ ariba run mydata Sample1_1.fq.gz Sample1_2.fq.gz output
```

## 结果分析




## Reference

1. [ariba](https://github.com/sanger/ariba)
