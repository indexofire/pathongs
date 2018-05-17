# Bioawk

[Bioawk](https://github.com/lh3/bioawk) 是 Heng Li 开发的 awk 扩展工具，增加了对压缩的 BED, GFF, SAM, VCF, FASTA/Q 等文件格式的支持，并内建一些函数，适用于NGS数据的快速输入输出。

- gc($seq) Returns the GC percentage of a sequence.
- meanqual($seq) Returns the average quality of the fastq sequence.
- reverse($seq) Returns the reverse of the sequence.
- revcomp($seq) Returns the reverse complement of the sequence.
- qualcount($qual, threshold) Returns the number of quality values above the threshold parameter.
- trimq(qual, beg, end, param) Trims the quality string qual in the Sanger scale using Richard Mott's algorithm (used in Phred). The 0-based beginning and ending positions are written back to beg and end, respectively. The last argument param is the single parameter used in the algorithm, which is optional and defaults 0.05.
- and(x, y) bit AND operation (& in C)
- or(x, y) bit OR operation (| in C)
- xor(x, y) bit XOR operation (^ in C)



## 1. 安装

```bash
$ sudo apt-get install bison
$ git clone https://github.com/lh3/bioawk
$ cd bioawk && make
$ sudo cp bioawk /usr/local/sbin
```

## 2. 使用

构建测序数据分析的 workflow 时，当 fastq 数据在做完 trimming 后，我们往往要关注剩下多少 reads，可以用 Bioawk 进行快速统计。

```bash
# 快速统计fastq里的reads数量
$ bioawk -c fastx 'END { print NR }' my_fastq.tar.gz
```

在一些特殊场合里，需要分析 reads 系列，用 Bioawk 可以很方便快速来完成。比如统计特殊碱基开头的 reads 数。

```bash
# 统计以GATTAC开头的reads的数量
$ bioawk -c fastx '$seq~/^GATTAC/ {++n} END { print n }' my_fastq.tar.gz
```

如果想看GC含量大于60%的reads的质量情况，可以用下面方法：

```bash
$ bioawk -c fastx '{if (gc($seq)>0.6) printf "%s\n%s\n\n", $seq, $qual}' my_fastq.tar.gz
```

生成reads的平均Q值

```bash
$  bioawk -c fastx '{print $name, meanqual($seq)}' my_fastq.tar.gz > meanqual.txt
```