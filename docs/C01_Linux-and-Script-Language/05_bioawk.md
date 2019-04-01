# bioawk

[bioawk](https://github.com/lh3/bioawk) 是 Heng Li 开发的 awk 扩展工具，增加了对压缩的 BED, GFF, SAM, VCF, FASTA/Q 等文件格式的支持，并内建一些函数，适用于NGS数据的快速输入输出。因为其是 awk 的功能扩展，所以既可以独立以 bioawk 命令，也可以将其 link 到 awk/gawk ，直接代替 awk。

### 内建函数

bioawk 内建了一些 awk 没有的专用于生物序列的函数：

- gc($seq) 获得序列的GC含量
- meanqual($seq) 获得序列的平均Q值
- reverse($seq) 将序列5'<->3'反向
- revcomp($seq) 将序列反向互补
- qualcount($qual, threshold) 获得大于阈值的序列中碱基数
- trimq(qual, begin, end, param) Trims the quality string qual in the Sanger scale using Richard Mott's algorithm (used in Phred). The 0-based beginning and ending positions are written back to beg and end, respectively. The last argument param is the single parameter used in the algorithm, which is optional and defaults 0.05.
- and(x, y) bit AND operation (& in C)
- or(x, y) bit OR operation (| in C)
- xor(x, y) bit XOR operation (^ in C)

### 内建变量

- bed: 1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts
- sam: 1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual
- vcf: 1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info
- gff: 1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute
- fastx: 1:name 2:seq 3:qual

## 1. 安装

在 Linux 中安装 bioawk

```bash
# Ubuntu 中源码编译安装
$ sudo apt-get install bison
$ git clone https://github.com/lh3/bioawk
$ cd bioawk && make
$ sudo cp bioawk /usr/local/sbin

# Archlinux aur 安装包安装
$ git clone https://aur.archlinux.org/bioawk.git
$ cd bioawk
$ makepkg -si

# conda 安装
(base)$ conda install bioawk
```

查看说明文档

```bash
$ cd bioawk
$ man ./awk.1
```

## 2. 使用示例

使用示例部分使用 NCBI SRA 数据库的 illumina 测序数据 SRR1175124 为例

- fastq: SRR1175124
- bed:
- vcf:

### 处理 fastq

**计算reads数量**

构建测序数据分析的 workflow 时，当 fastq 数据筛选后，如果要查看剩下多少 reads，可以用 bioawk 进行快速统计。

```bash
# 去除低质量reads后，快速统计fastq里的reads数量
$ bioawk -c fastx 'END { print NR }' SRR1175124_1.fastq.gz
102354
$ echo `zcat SRR1175124_1.fastq.gz | wc -l` / 4 | bc
102354

# 统计以GATTAC开头的reads的数量
$ bioawk -c fastx '$seq ~ /^GATTAC/ {++n} END { print n }' SRR1175124_1.fastq.gz
# 用 grep 实现
$ zcat SRR1175124_1.fastq.gz | grep '^GATTAC' | wc -l

# 查看GC含量大于60%的reads的质量情况
$ bioawk -c fastx '{if (gc($seq)>0.6) printf "%s\n%s\n\n", $seq, $qual}' \
> SRR1175124_1.fastq.gz

# 生成reads的平均Q值
$ bioawk -c fastx '{print $name, meanqual($seq)}' \
> SRR1175124_1.fastq.gz > meanqual.txt

# 将前5条reads序列的5'->3'变为3'->5'
$ bioawk -c fastx '{if(NR<6) print reverse($seq)}' SRR1175124_1.fastq.gz

# 输出前5条reads序列的反义互补序列
$ bioawk -c fastx '{if(NR<6) print revcomp($seq)}' SRR1175124_1.fastq.gz

# 输出第100条reads的信息
$ bioawk -c fastx 'NR==100{print ">"$name, gc($seq)}' SRR1175124_1.fastq.gz

# 统计SRR1175124_1.fastq.gz中reads碱基Q值高于35的个人小于50的reads数
$ bioawk -c fastx '{if(qualcount($qual, 35)<50) print $name,$seq}' \
> SRR1175124_1.fastq.gz | wc -l

# 双末端测序的碱基Q值高于35的个数小于50的reads
$ bioawk -c fastx '{if(qualcount($qual, 35)<50) print $name}' \
> SRR1175124_1.fastq.gz SRR1175124_2.fastq.gz | sort | uniq -c | \
> awk '$1==2{print}'

# 输出双末端测序的碱基Q值高于35的个数小于50的reads中以5’-ATG的reads
$ for i in $(bioawk -c fastx '{if(qualcount($qual, 35)<50) print $name}' \
> SRR1175124_1.fastq.gz  SRR1175124_2.fastq.gz | sort | uniq -c | \
> awk '$1==2{print $2}'); do bioawk -c fastx -v i="$i" \
> '{if($name==i && $seq ~ /^ATG/) print $name"\n"$seq"\n"$qual"\n"}' \
> SRR1175124_1.fastq.gz SRR1175124_2.fastq.gz; done
```

```bash
$ bioawk -c bed '{ print $end - $start }' .bed
```

### 处理 bam

### 处理 vcf

### 处理bed

### 处理gff


