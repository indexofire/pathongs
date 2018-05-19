# 炭疽芽胞杆菌基因组 SNPs 分析

疾控部门属于保障单位，需要加强对生物恐怖事件的防控和应对处置能力。最著名的细菌性生物恐怖战剂就是含有炭疽芽胞杆菌的白色粉末。日常工作中无法接触到这些病原，因此在经验和技能方面较为缺乏。在G20期间，内部也进行了相关项目的应急拉练，除此之外为了提升能力建设，应对可能发生的疑似事件，在实验室快速分析和溯源工作做好能力储备，我们对炭疽芽胞杆菌的基因组数据进行了基本分析，对相关资料和信息进行了储备。

主要通过构建相应的分析流程，获得菌株基因组数据后可以第一时间就进行溯源方面的研究。从而在这准备工作的过程中形成了此文档。主要目的也是为了抛砖引玉，做到举一反三，给其他细菌性病原微生物基因组溯源分析时提供参考。

软件中使用的一些参数是根据自己的服务器设置的，具体运行时需要根据计算机配置进行相应的调整。我们服务器的配置如下：

!!! note
    * CPU: E3 2650v3 x2
    * Memory: ECC 16G x8
    * Harddisk: 7200转 1T x4 组 Raid 5 + 1 Hotspare
    * System: Ubuntu 16.04 LTS 64bit

## 1. 获取公共数据库数据

数据采用 NCBI SRA 数据库中的公共数据来作分析。首先检索 NCBI SRA 数据库中所有使用 illumina 平台对炭疽芽胞杆菌基因组进行 PE 测序的 大于 100M 的 NGS 数据，并在下载这些 sra 数据后转换成 fastq.gz 格式。使用的工具主要为 `edirect` 和 `sra-tools`。

```bash
$ esearch -db sra -query '"Bacillus anthracis"[Organism] AND \
> (Hiseq[ALL] OR Miseq[ALL]) AND "strategy wgs"[Properties]' | \
> efetch -db sra -format runinfo | grep 'PAIRED' | grep '^[DSE]RR*' | \
> awk -F',' '{if($5>100000000) print $1}' | sort -n | \
> xargs fastq-dump --split-files --gzip --outdir bacillus_anthracis
```

虽然我们检索时限定了 PE 测序数据，但是很多时候数据库里的信息不准确，还是会有很多 SE 的数据会被下载。确定 PE 或者 SE 的方法很多。

```bash
# sra-stat 统计 sra 文件，nreads=1 表示 SE 测序，nreads=2 表示 PE 测序。程序统计整个文件，速度比较慢
$ sra-stat -xs *.sra | grep "nreads"
# 批量删除
$ parallel 'sra-stat -xs' ::: *.sra | xtract -Pattern Run -element Run@accession \
> Statistics@nreads | awk '{if($2==1) print $1}' | xargs rm

# fastq-dump 输出一个 reads，计算行数如果是4，那么就是 SE，如果是8,则是 PE
$ fastq-dump -X 1 --split-spot -Z SRR955386.sra | wc -l

# fastq-dump 自己判断，当 sra 文件是 SE 测序时，fastq-dump 只能生成1个 *_1.fastq 文件
$ fastq-dump --split-files ERR493452.sra

# fastq-dump 只能调用单核转换，可以使用 parallel 工具并行处理加快速度。
$ parallel "fastq-dump --split-files --gzip --outdir bacillus_anthracis" ::: *.sra
```

有时候有些sra文件在 dump 时会出错，不能生成对应的 PE fastq 文件，这就需要比较后删除了。如果是 mate pair 的数据，会生成 `*_3.fastq.gz`，如果不需要则可手工删除。

```bash
# 核对 PE 双端测序文件是否一致
$ ls -l *_1.fastq.gz | awk '{print $9}' | awk -F'_' '{print $1}' > R1.txt
$ ls -l *_2.fastq.gz | awk '{print $9}' | awk -F'_' '{print $1}' > R2.txt
$ diff R1.txt R2.txt
# diff 如果2个文件内容一致，则终端不会有输出。否则就要删除文件
```

## 2. 数据的前期处理

由于是公共数据库下载的数据，并不能保证测序实验质量或者数据提交者是否提交的是质控后的数据。因此建议要对数据做一些前期的质控处理。

### 2.1 数据 QC

首先对基因组 GC 含量和 Q 值进行初步筛选，对于偏差较大的数据考虑直接剔除（或者用其他软件验证看是否是错误物种）。这里使用的工具为 `bioawk` 或 `parallel`。

```bash
# 因为 reads 覆盖度可能并不均一，高覆盖度的区域 GC 含量比重会略高。但是大部分情况下，正常测序的基因组这种区域相对来说不多，平均到基因组后会整体 GC 含量影响不大。
$ for i in *.fastq.gz; do bioawk -c fastx 'BEGIN{n=0;q=0}{n+=gc($seq);q+=meanqual($seq)}END{print $name,n/NR,q/NR}' $i >> result.txt; done
# awk 类工具是单进程的，为了加速可以使用 parallel 来并行计算
$ parallel "bioawk -c fastx 'BEGIN{n=0;q=0}{n+=gc(\$seq);q+=meanqual(\$seq)}END{print \$name,n/NR,q/NR}' >> result.txt" ::: *.fastq.gz
# 按照 GC 含量排序
$ cat result.txt | sort -n -k3
```

R 绘制 gc 分布图，如果 gc 含量偏差超过 10% 时就剔除该基因组数据。

```r
> library(ggplot2)
> data<-read.table("result.txt")
> qplot(v2, v3, data=data)
```

结果显示有一个黑点的GC含量特别低，且Q值特别高。结果显示测序实验 SRR2155551 和 SRR2164197 的 GC 含量仅为22%，而其他样本 GC 含量范围在34%～41%。且这2个实验的样本均为同一个，是1个样本的2此测序（或者数据上传了2次），因此将其剔除。

### 2.2 去除接头

其次要去除接头的污染。因为分析流程中不仅包括 mapping 的方式，还包含 de novo assembly，为了避免接头序列对基因组拼接的影响，这里最好进行。这里使用的工具是 `FastQC` 和 `fadapa`

```
# 使用脚本 scan_adaptors.py 来扫描下载的高通量基因组测序数据是否有接头污染的情况。
$ fastqc -t 40 -q --extract *.fastq.gz
$ python scan_fastqc_report.py -d qc
```

```python
#!usr/bin/env python
# -*- coding: utf-8 -*-
# scripts name: scan_fastqc_report.py
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from fadapa import Fadapa
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', dest='qcfile', type=str, help="fastqc extract data file")
parser.add_argument('-o', '--output', action='store', dest='adaptor', type=str, help="adaptors file")
args = parser.parse_args()
qc_file = Fadapa(args.qcfile)

if qc_file.summary()[-3][0] == 'fail':
    adaptors = []
    ors = qc_file.clean_data("Overrepresented sequences")
    ors.pop(0)

    for (index, seq) in enumerate(ors):
        if 'NNNNNN' in seq[0]:
            continue
        adaptors.append(SeqRecord(Seq(seq[0], generic_nucleotide), id="adaptor_%d" % (index+1)))

    if not os.path.exists(args.adaptor):
        os.makedirs(args.adaptor)
    SeqIO.write(adaptors, "%s/adaptor.fasta" % args.adaptor, "fasta")
    print "Overrepresented sequences has been save to adaptors.fasta"
else:
    print "No Overrepresented sequences"
```


## 3. Pangenome 分析

### 3.1 基因组组装

```
$ spades.py
```

### 3.2 基因组注释

```
$ prokka
```

### 3.3 Roary

```
$ roary -p -e mafft *.gff
```

### 3.4 Harvest

```
$
```

## 4. SNPs 构建进化树

### 4.1 基于参考基因组 Mapping 的方法

**Snippy**

```
# snippy mapping
$ for i in *.fastq.gz | sort | uniq; \
> do snippy --cpus 20 --outdir $i -ref reference.fa \
> --R1 $i*_1*.fastq.gz --R2 $i*_2*.fastq.gz; done

# snippy-core 汇集 snps
$ snippy-core --prefix core-snps SRR* DRR* ERR*

# raxml 绘制进化树
$ raxml -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s core-snps.aln -n ex -T 40

# 将 raxml 进化树结果文件从服务器传到本地电脑，用 figtree 生成图片
$ scp user@server-ip:/path/RAxML_bestTree.ex .
$ figtree RAxML_bestTree.ex &
```

**Tools**


## 5. 基于组装结果的核心基因组序列


## 相关软件安装

**edirect**

```bash
$ wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
$ sudo tar zxf edirect.tar.gz -C /opt
$ sudo chown -R root:root /opt/edirect
$ cd /opt/edirect && sudo ./setup.sh
$ sudo ln -s /opt/edirect/* /usr/local/sbin/
```

**sra-tools**

```bash
$ wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-ubuntu64.tar.gz
$ sudo tar zxf sratoolkit.2.7.0-ubuntu64.tar.gz -C /opt/sratoolkit
$ sudo chown -R root:root /opt/sratoolkit && cd /opt/sratoolkit
$ sudo ln -s `pwd`/bin/* /usr/local/sbin/
```

**bio-awk**

```bash
$ sudo apt-get install bison
$ git clone https://github.com/lh3/bioawk
$ cd bioawk && make
$ sudo cp bioawk /usr/local/sbin
```

**parallel**

```bash
$ sudo apt install parallel
```

**gplot**

```bash
$ sudo apt install gnuplot
$ git clone https://github.com/RhysU/gplot
$ sudo cp gplot/gplot /usr/local/sbin
```

**FastQC**

```bash
$ wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
$ sudo unzip fastqc_v0.11.5.zip
$ sudo mv fastqc_v0.11.5 /opt/fastqc && sudo chown -R root:root /opt/fastqc
$ sudo cp /opt/fastqc/fastqc /usr/local/sbin
```

**Fadapa**

```bash
$ sudo pip install fadapa
```

**snippy**

```bash
$ git clone https://github.com/tseemann/snippy
$ sudo cp snippy/bin/* /usr/local/sbin
$ sudo cp snippy/binaries/linux/* /usr/local/sbin
```

**RAxML**

```bash
$ wget https://github.com/stamatak/standard-RAxML/archive/v8.1.17.tar.gz
$ sudo tar zxf v8.1.17.tar.gz -C /opt/raxml
$ sudo chown -R root:root /opt/raxml && cd /opt/raxml
$ sudo make -f Makefile.SSE3.PTHREADS.gcc
$ sudo ln -s `pwd`/raxmlHPC-PTHREADS-SSE3 /usr/local/raxml
```

## Reference

1. edirect help book: http://www.ncbi.nlm.nih.gov/books/NBK179288/
2. sra-tools wiki: https://github.com/ncbi/sra-tools/wiki
3. parallel tutorial: https://www.biostars.org/p/63816/
4. gnuplot tutorial: http://lzz5235.github.io/2016/01/12/gnuplot.html
5. gplot usage: https://github.com/RhysU/gplot/blob/master/README.md
