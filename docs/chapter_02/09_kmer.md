# k-mer 分析

[khmer][]

```bash
$ pip install khmer
```

### 通过 conda 安装

用 conda 新建一个虚拟环境安装 [khmer][]。当前版本为 [2.12][http://khmer.readthedocs.io/en/v2.1.2/]，安装过程还会自动 [screed][https://github.com/dib-lab/screed]（一个轻量级 python 序列数据库工具）。

```bash
$ conda create -n khmer
$ conda activate khmer
(khmer)$ conda install khmer ipython
```

[khmer][] 有一套命令行工具，可以非常方便的运行使用。另外也可以使用 [khmer][] 的 API，编写 python 代码调用。

**load-into-counting.py**

```bash
$ load-into-counting.py -T 4 -k 32 -x 1e8 -N 4 counts R1.fastq.gz
```


**normalize-by-median.py**

去除冗余序列。脚本程序将 Fast[A/Q] 格式数据文件列表化，

脚本输出的序列可以被组装软件比如 SPAdes 等进行序列组装。

```bash
$ normalize-by-median.py -p -k 21 -o - read.fq output.fq
```

**readstats.py**

统计 FASTA/Q 序列信息。

```bash
$ readstats.py -o R1.output R1.fastq.gz
$ cat R1.output

found 449472737 bps / 1879073 seqs; 239.2 average length -- R1.fastq.gz
---------------
449472737 bp / 1879073 seqs; 239.2 average length -- R1.fastq.gz
---------------
449472737 bp / 1879073 seqs; 239.2 average length -- total

```


## khmer 用途

对 de novo assembly 进行数字均一化进行预处理

## Python API

```bash
$ ipython
In [1]: import khmer
In [2]: counts = khmer.Counttable(31, 1e7, 4)
In [3]: counts.consume_seqfile('R1.fastq.gz')
(25000, 1150000)
# 31 是 kmer 长度，1e7 是 table size，4 是 table 数
# 25000 表示序列数量，高通量测序里即 reads 数。
# 1150000 表示指定 kmer 切割后特定 kmer 序列数
In [4]: counts.get('TGACTTTCTTCGCTTCCTGACGGCTTATGCC')
105
```


## Reference
1. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4608353/
2. http://khmer.readthedocs.io/en/latest/

[khmer]: https://github.com/dib-lab/khmer/ "KhmerM"
