# Minimap & miniasm


## 1. 安装

```bash
# 在之前建立的 nanopore 虚拟环境中安装
(nanopore)$ conda install minimap
(nanopore)$ conda install miniasm
```

## 2. 使用

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
(seq)$ edirect -db 'sra' -query 'nanopore salmonella' | efetch -db

(seq)$ prefetch -v

(seq)$ fastq-dump

(nanopore)$ minimap -Sw5 -L100 -m0 -t40
(nanopore)$ miniasm -f reads.fq
```
