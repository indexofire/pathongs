# 长读长测序数据组装

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

目前 Pacbio 和 Nanopore 测序技术可以成生 K 级别甚至 M 级别的reads，但是其目前来说 base calling 产生的错误还是较高，因此一般称其 long noisy reads 数据。本节介绍 pacbio/nanopore 测序数据进行基因组组装。用到的软件有：

- minimap2: all vs all 比对
- miniasm: 组装
- racon: 矫正

minimap 是 Li heng 针对 long noisy reads 的比对工具，目前已经更新的了 v2 版本。miniasm 也是其开发的用于组装的工具。

## 1. 搭建分析环境

```bash
# 在之前建立的 nanopore 虚拟环境中安装
$ conda create -n nanopore minimap miniasm
```

下载实验室数据

```bash
# SRR8494940 是一个使用 minION 测序的大肠埃希菌 ST131-H22 菌株
(seqs)$ prefetch -v SRR8285369
(seqs)$ fastq-dump --gzip SRR8285369.sra
```

开始分析之前的数据文件

```
.
├── SRR8285369.fastq.gz
└── SRR8285369.sra
```

```bash
# 使用 minimap 进行 all vs all 自我比对
(nanopore)$ minimap -Sw5 -L100 -m0 -t 4 SRR8285369.fastq.gz SRR8285369.fastq.gz | gzip -1 > SRR8285369.paf.gz
# 使用 minimap2 进行 all vs all 自我比对
(nanopore)$ minimap2 -x ava-ont -t 4 SRR8285369.fastq.gz SRR8285369.fastq.gz | gzip -1 > SRR8285369.paf.gz

# miniasm 组装基因组
(nanopore)$ miniasm -f SRR8285369.fastq.gz SRR8285369.paf.gz > SRR8285369.gfa
# 生成 fasta 个是序列
(nanopore)$ awk '$1 ~ /S/ {print ">"$2"\n"$3}' SRR8285369.gfa > SRR8285369.fasta

# 使用 racon 矫正
# round 1
(nanopore)$ minimap2 -t 4 SRR8285369.fasta SRR8285369.fastq.gz > correct1.gfa.paf
(nanopore)$ racon -t 4 SRR8285369.fastq.gz correct1.gtf.paf SRR8285369.fasta > racon1.fasta
# round 2 （可选）
(nanopore)$ minimap2 -t 4 racon1.fasta SRR8285369.fastq.gz > correct2.gtf.paf
(nanopore)$ racon -t 4 SRR8285369.fastq.gz correct2.gtf.paf racon1.fasta > racon2.fasta

(nanopore)$ mv racon2.fasta SRR8285369_corr.fasta
(nanopore)$ mash SRR8285369.fasta SRR8285369_corr.fasta
```

组装完成后的数据文件

```
.
├── correct1.gfa.paf
├── correct2.gfa.paf
├── racon1.fasta
├── racon2.fasta
├── SRR8285369_corr.fasta
├── SRR8285369.fasta
├── SRR8285369.fastq.gz
├── SRR8285369.gfa
├── SRR8285369.paf.gz
└── SRR8285369.sra
```

`long-noisy-assembly.sh` 创建一个 shell 脚本，便于组装操作:

```bash
#!/bin/bash
source activate nanopore

cpus = `cat /proc/cpuinfo | grep "processor" | wc -l`
name = `awk -F'.' '{print $1}'`

minimap2 -x $1 -t $cpus $2 $2 | gzip -1 > ${name}.paf.gz
miniasm -f $2 ${name}.paf.gz > ${name}.gfa
awk '$1 ~ /S/ {print ">"$2"\n"$3}' ${name}.gft > ${name}.fasta

minimap2 -t $cpus ${name}.fasta $2 > ${name}_corr1.gfa.paf
racon -t $cpus $2 ${name}_corr1.gfa.paf ${name}.fasta > ${name}_racon1.fasta
minimap2 -t $cpus ${name}_racon1.fasta $2 > ${name}_corr2.gfa.paf
racon -t $cpus $2 ${name}_corr2.gfa.paf ${name}_corr2.gtf.paf > ${name}_racon2.fasta
cp ${name}_racon2.fasta ${name}_corr.fasta
```

## 2. 使用

### 2.1 miniasm 使用

miniasm: https://github.com/lh3/miniasm

miniasm 是一个基于 Overlap OLC 算法的 de novo assembler，可以对类似 pacbio 或 nanopore 技术获得的长高错误率的序列进行组装。miniasm 采用自我 reads 的 all vs all 比对策略，

```bash
# 组装 pacbio reads 数据
(nanopore)$ miniasm -x ava-pb -t 4

# 组装 nanopore 数据
(nanopore)$ miniasm -x ava-ont -t 4
```

miniasm 主要参数

- `R`
- `m` 最小匹配序列长度，默认值100
- `i`

### 2.2 minimap 使用

```bash
# 多对多比对
(nanopore)$ minimap -Sw5 -L100 -m0 -t40 reads.fq reads.fq | gzip -1 > reads.paf.gz

# 组装
(nanopore)$ miniasm -f reads.fq reads.paf.gz > reads.gfa

# 将 gfa 格式转换为 fasta 格式
(nanopore)$ awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > reads.fa

# 拼接结果的 contigs 数
(nanopore)$ grep ">" reads.fa | wc -l
```

## 3. 实例

```bash
# 使用 edirect 获取 NCBI sra 数据库中的 nanopore 测序的沙门菌
(seq)$ edirect -db 'sra' -query 'nanopore salmonella' | efetch -db

# 下载数据
(seq)$ prefetch -v

# 将 sra 格式的数据转换成 fastq 格式
(seq)$ fastq-dump

# 生成模拟数据
(nanopore)$ minimap -Sw5 -L100 -m0 -t40
(nanopore)$ miniasm -f reads.fq
```

## Reference

https://yiweiniu.github.io/blog/2018/03/Genome-assembly-pipeline-miniasm-Racon/
https://github.com/arvindsundaram/IN-BIOSx000/blob/2018/Assembly/practicals/07_Assembly_using_minasm+racon.md
