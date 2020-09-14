# IGV

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    igv是最常用的比对结果可视化软件之一，使用java开发可以跨平台使用。

## 1. 准备数据

```bash
# 下载参考基因组 NC_012967
$ esearch -db nuccore -query 'NC_012967' | efetch -format fasta > NC_012967.1.fasta

# 下载比对序列sra
$ prefetch -v SRR030257
```

## 2. 安装

```bash
$ wget -O http://iubio.bio.indiana.edu/soft/molbio/readseq/java/readseq.jar
$ java -cp ~/app/readseq/readseq.jar run NC_012967.1.gbk -f GFF -o NC_012967.1.gff
$ wget https://github.com/broadinstitute/IGV/archive/v2.3.37.zip
$ unzip v2.3.37.zip
$ cd v2.3.37
$ java -Xmx2g -jar igv.jar
```

!!! infomation
    Linux下使用igv建议采用java8，其他版本可能会出错。因此通过conda安装也是一个不错的方式。

## 3. 生成比对文件

首先将 sra 格式的比对文件转换成 fastq 格式，并进行质控去除接头和外缘污染等序列。然后用samtools

**用 samtools 转换数据格式**

```bash
$ samtools faidx NC_012967.1.fasta
$ samtools view -b -S -o bowtie/SRR030257.bam bowtie/SRR030257.sam
$ samtools sort bowtie/SRR030257.bam bowtie/SRR030257.sorted
$ samtools index bowtie/SRR030257.sorted.bam
```

**分析结果**

弹出的IGV图形界面，下拉菜单选择`Genomes`/`Create .genome File...`,在弹出的窗口中`FASTA file`选项选择基因组fasta文件，在`Gene File`选项选择转换的文件。

## 参考资料
