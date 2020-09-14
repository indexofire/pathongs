# Seqkit

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

[Seqkit]() 是由国人 [shenwei]() 开发，功能强大的命令行序列处理工具套件。

## 使用

| 子命令 | 用途 |
| ------ | ---- |
| common | find common sequences of multiple files by id/name/sequence |
| concat | concatenate sequences with same ID from multiple files |
| convert | convert FASTQ quality encoding between Sanger, Solexa and Illumina |
| duplicate | duplicate sequences N times |
| faidx | create FASTA index file and extract subsequence |
| fq2fa | convert FASTQ to FASTA |
| fx2tab | convert FASTA/Q to tabular format (with length/GC content/GC skew) |
| genautocomplete | generate shell autocompletion script |
| grep | search sequences by ID/name/sequence/sequence motifs, mismatch allowed |
| head | print first N FASTA/Q records |
| help | Help about any command |
| locate | locate subsequences/motifs, mismatch allowed |
| mutate | edit sequence (point mutation, insertion, deletion) |
| range |  print FASTA/Q records in a range (start:end) |
| rename | rename duplicated IDs |
| replace | replace name/sequence by regular expression |
| restart | reset start position for circular genome |
| rmdup |  remove duplicated sequences by id/name/sequence |
| sample | sample sequences by number or proportion |
| seq | transform sequences (revserse, complement, extract ID...) |
| shuffle | shuffle sequences |
| sliding | sliding sequences, circular genome supported |
| sort | sort sequences by id/name/sequence/length |
| split | split sequences into files by id/seq region/size/parts (mainly for FASTA) |
| split2 | split sequences into files by size/parts (FASTA, PE/SE FASTQ) |
| stats |  simple statistics of FASTA/Q files |
| subseq | get subsequences by region/gtf/bed, including flanking sequences |
| tab2fx | convert tabular format to FASTA/Q format |
| translate |  translate DNA/RNA to protein sequence |
| version |  print version information and check for update |



### grep 命令

spades 组装完毕的序列，提取其中某个 node 或者某些 nodes 序列：

```bash
$ seqkit grep -r -p ^node_1$ scaffolds.fasta > node_1.fasta
```

将 scaffolds 中所有 nodes 提取到单独的 fasta 文件中：

```bash
$ for i in $(awk '/>/' scaffolds.fasta | awk '{print $1}'); \
> do seqkit grep -p ${i:1} scaffolds.fasta > ${i:1}.fasta; \
> done
```

### seq 命令

```bash
$ seqkit seq -h

```

```bash tab="-n和-i"
# 将 reads 的名称信息输出，如果要和 bioawk $name 相同，使用 -i 参数
$ seqkit seq -n SRR1175124_1.fastq.gz SRR1175124_2.fastq.gz
SRR1175124.102354 102354 length=150
...
$ bioawk -c fastx '{print $name}' SRR1175124_1.fastq.gz SRR1175124_2.fastq.gz
SRR1175124.102354
...
$ seqkit seq -in SRR1175124_1.fastq.gz SRR1175124_2.fastq.gz
SRR1175124.102354
...

#
```

```bash tab="-l和-u"
# 输出小写碱基
$ seqkit seq -sl SRR1175124_1.fastq.gz
# 输出大写碱基
$ seqkit seq -su SRR1175124_1.fastq.gz
```




## 使用案例

```bash
# 下载 NC_001477
$ efetch -db nuccore -id NC_001477 -format fasta > NC_001477.fasta
# 翻译 DNA 序列为氨基酸序列
$ seqkit translate -j 4 -o NC_001477.pep NC_001477.fasta
#
```

!!! note "对注释的多拷贝基因提取序列后进行序列比对"

```bash
# 用prokka注释基因组拼接序列
$ prokka --rfam --prefix sample --outdir sample contig.fa
# 用 seqkit 提取多拷贝基因 ncRNA Qrr 的序列
$ grep 'Qrr' sample/sample.gff | for i in $(awk '{print $9}'); do seqkit grep -p ${i:3:14} sample.ffn >> qrr.fa; done
# mafft 序列比对
$ mafft qrr.fa > qrr.mafft
# raxml 快速构建进化树
$ raxmlHPC -f a -p 12345 -x 12345 -# 1000 -m GTRGAMMA -s qrr.mafft -n qrr
```
