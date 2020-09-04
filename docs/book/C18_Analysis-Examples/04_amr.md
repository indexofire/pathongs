# 基于高通量测序技术进行病原细菌耐药研究

---

!!! Abstract "内容简介"
    本节介绍基因组数据进行耐药基因的研究实践。

## 数据质控

```bash
# trimmomatic 取出低质量序列和接头序列
# 去除 NexteraPE-PE.fa 文件中的接头文件
# 取出头尾部3碱基范围中Q值低于3的质量碱基
# 以4个碱基长度为读码长度，取出平均质量低于15的剪辑
# 取出长度小于36bp的碱基
$ export ADAPTER=$CONDA_PREFIX/share//trimmomatic-0.39-1/adapters
$ trimmomatic PE -threads 4 -trimlog 02_trim/trim.log -summary 02_trim/ \
> summary.log DRR015013_1.fastq.gz DRR015013_2.fastq.gz 02_trim/ \
> DRR015013_1.clean.fq.gz 02_trim/DRR015013_1.trim.fq.gz 02_trim/ \
> DRR015013_2.clean.fq.gz 02_trim/DRR015013_2.trim.fq.gz \
> ILLUMINACLIP:$ADAPTER/NexteraPE-PE.fa:2:30:10 \
> LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

$ fastp
```

## 基因组数据处理

### 1.基因组组装

### 2.基因组比对

### 3.基因组注释

swiss-prot
KEGG
COG
GO
ncRNA

### 4.基因组特征展示

circos

### 5.细菌


## 基因组比较

### 1.共线性分析


## 多个基因组分析
