# 用 awk 处理 NGS 数据

---

![awk](https://pbs.twimg.com/media/DeLcVfSWAAAw6OZ.jpg)

!!! Abstract "内容简介"
    awk 工具常用来对文本文件的内容按行进行匹配搜索，并对搜索条件匹配的内容进行操作。高通量测序技术带来大量的数据，以及初级分析生成的各种按列的数据结果，awk 就可以大显身手为测序数据的转换与分析提供方便的操作。

---

**我们首先来看一个最简单的 awk 命令程序：**

```bash
# 序列文件contigs.fasta文件包含2条DNA序列，分别是seq1和seq2。
$ cat contigs.fasta
>seq1
ATGGGGTCCGATTACCG
>seq2
CGGNTTACCGYGACAGT

# 提取2条序列的id
$ awk '/>/' contigs.fasta
>seq1
>seq2

# 提取ATG开始的序列
$ awk '/^ATG/' contigs.fasta
ATGGGGTCCGATTACCG
```

!!! note
    awk 用模式匹配时，与 grep 类似。使用这个命令很简单，就是将 fasta 文件中包含 '>' 字符的行输出，实际上就是将所有的序列头内容输出。作用类似`grep ">"`。

**利用 awk 可以很方便的调整或输出列：**

```bash
# 例如这样数据文件input_data
$ cat input_data
gene1 100-300 description
gene2 1000-1500 description
...

# 当第一列包含gene1字符时，打印第二列信息基因在基因组上的位置
$ awk '$1 ~/gene1/ {print $2}' input_data
100-300

# 用 grep 和 sed 实现类似效果
$ cat input_data | grep gene1 | sed {$1}

# 也可以很方便的调整列的顺序
$ awk '{print $3"\t"$1"\t"$2}' input_data
description gene1 100-300
description gene2 1000-1500
...
```

---

## 1. awk 命令构成

awk 程序主体主要包括 pattern，action和input file：

* pattern 表示所要搜索的内容，可以用正则表达式。
* { action } 则表示搜索匹配后要做的操作。
* input file：所要搜索的输入内容，比如前面例子的 contigs.fasta

!!! info "知识点"
    awk可以在shell里直接运行：

    ```bash
    $ awk 'pattern { action }' input_file
    ```

    也可以写成awk程序来运行。一个 awk 程序往往如下所示：

    ```awk
    #!/usr/bin/awk -f
    pattern { action }
    ```

awk 可以不需要输入文件；对于pattern和action来说，2者至少要有一个才能运行。如果没有pattern，则默认匹配任何输入，按行输出并执行action。如果没有action，则匹配pattern并按行输出不做额外操作。

```bash
# 尝试以下语句，在终端打印`hello world`。
$ awk 'BEGIN { print "hello world" }'
# 这里BEGIN是pattern，{ print "hello world" } 是action

# 模拟cat输出终端输入的字符。
$ awk '{ print }'
# 这里省略了pattern，只有{print}做为action
```

上面这个例子可以保存成文件（如`hello_world.awk`）来运行，文件内容如下：

```awk
#!/usr/bin/awk -f
# 注意不同发行版的linux，awk路径有所区别

BEGIN { print "hello world" }
```

运行这个程序：

```bash
$ chmod +x hello_world.awk
$ ./hello_world.awk
```

## 2. awk 基本语法

### 2.1 Pattern

#### 2.1.1 正则表达式

pattern 可以采用正则表达式来匹配内容。用`/ regular expression /`来表示。

几个特殊的Patterns：

- BEGIN: awk开始输入内容之前，可以用来定义初始变量等
- END: awk结束输入内容之后，可以用来输出统计信息等
- BEGINFILE: awk输入内容为多个文件时，统一定义多个文件的BEGIN块
- ENDFILE: awk输入内容为多个文件时，统一定义多个文件的END块

#### 2.1.2 语句表达

```bash
# 当第一列内容是gene时，输出第二列内容。可用于注释文件
$ awk '$1 == "gene" { print $2 }' genes.list
```

### 2.2 Action

#### 2.2.1 if...else

awk的action支持条件判断语句`if else`

```bash
# if 语句的基本用法
$ awk '{if(condition) print $1}' input

# if else block
$ awk '{if(condition)print $1;else print $2}' input
```

#### 2.2.2 while

awk的action支持循环语句，可以用while控制。while 控制语句要换行，用4个空格划分控制块。

```bash
# 将会将列内容块只输出前3列内容，且没一个数据占据一行
$ awk '{i=1;while(i<=3){print $i;i++}}' input
```

### 2.3. 内建变量

#### 2.3.1 控制 awk 的变量

- FIELDWIDTHS：设置分割长度，如FS参数设置空格为分割符，可以用FIELDWIDTHS设置空格长度。
- FS：设置的分割符，默认是空格" "
- OFMT：默认格式化数字，默认是%.6g，
- OFS：设置print输出各字段的符号，默认是空格" "
- ORS：设置print输出多行的符号，模式是\n
- RS：设置行的分割符，默认为换行符\n

### 2.3.2 传递信息的变量

- ARGC：awk命令传递的参数个数
- ARGV：以数组形式保留awk传递的所有参数
- ENVIRON：用来在awk中调用环境变量
- FILENAME：当前正在处理的文件的名称
- FNR：当前文件awk运行的行号
- NF：表示当前行中有多少个字段
- NR：表示当前处理了多少行，当多个文件处理时，该值是FNR叠加的
- RLENGTH：表示匹配模式命中的匹配长度

### 2.4. 内建函数

#### 2.4.1 数学函数

可以直接调用awk的函数进行计算

- atan2(y, x)
- cos(x)
- exp(x)
- int(x)
- log(x)
- rand()
- sin(x)
- sqrt(x)
- srand([ x ])

#### 2.4.2 字符串函数

- asort(source [ , dest [ , how ] ] ) #
- asorti(source [ , dest [ , how ] ] ) #
- gensub(regexp, replacement, how [ , target ] ) #
- gsub(regexp, replacement [ , target ] )
- index(in, find)
- length( [ string ] )
- match(string, regexp [ , array ] )
- patsplit(string, array [ , fieldpat [ , seps ] ] ) #
- split(string, array [ , fieldsep [ , seps ] ] )
- sprintf(format, expression1, ...)
- strtonum(str) #
- sub(regexp, replacement [ , target ] )
- substr(string, start [ , length ] )
- tolower(string)
- toupper(string)

#### 2.4.3 输入输出函数
- close(filename [ , how ] )
- fflush( [ filename ] )
- system(command)

#### 2.4.4 时间函数
- mktime(datespec)
- strftime( [ format [ , timestamp [ , utc-flag ] ] ] )
- systime()

#### 2.4.5 Bit操作函数
- and(v1, v2 [, ...])
- compl(val)
- lshift(val, count)
- or(v1, v2 [, ...])
- rshift(val, count)
- xor(v1, v2 [, ...])

#### 2.4.6 其他
- isarray(x)
- bindtextdomain(directory [ , domain ] )
- dcgettext(string [ , domain [ , category ] ] )
- dcngettext(string1, string2, number [ , domain [ , category ] ] )

## 3. 示例

示例部分介绍 awk 在日常中，特别是处理 ngs 数据时的一些例子。

### 3.1 对列的操作

**调整列顺序**：在有GUI的操作系统里，一般采用类似 excel, calc 之类的软件导入数据文件，然后剪切各列调整顺序。如果用 awk 来解决也很方便，你只需要考虑好调整的各列顺序即可，action 里的{ print ...}顺序就是重新调整后的各列顺序：

```bash
$ awk '{ print $3, $5, $7, $2, $1, $4, $6 }' infile.txt > outfile.txt
```
**插入列**：有时候要插入一列数据，用awk可以很方便的实现：

```bash
$ awk '{ print $1, $2 "gene expression", $3}' infile.txt > outfile.txt
```

###  3.2 对行的操作

```bash
# 有时候数据里含有重复的行，而当你只需要唯一性数据时，就可以用这行程序，只保留具有唯一性的数据行。
$ awk '!x[$0]++' infile.txt > outfile.txt

# 不显示重复的行
$ awk '_[$0]++' infile.txt
# 当字段1重复时，打印整行
$ awk '_[$0]++{print $2}' infile.txt
```

### 3.3 格式化输出

```bash
# 格式化输出/etc/password内容，输出第1行到20行内容
$ awk -F: 'BEGIN{printf "%-8s %-3s", "User", "UID"}NR==1,NR==20{printf "%-8s %-3d", $1, $3}' /etc/passwd
```

### 3.4 合并文件

```bash
# 搜索当前目录下所有包含序列ATCCGGA的fasta文件，并输出文件名和行号以及序列
$ find . -type f -name "*.fasta" -exec awk '/ATCCCGGA/ {print FILENAME, NR, $0}' {} \+
```

---

## bioawk

[bioawk](https://github.com/lh3/bioawk) 是 Heng Li 开发的 awk 扩展工具，增加了对压缩的 BED, GFF, SAM, VCF, FASTA/Q 等文件格式的支持，并内建一些函数，适用于NGS数据的快速输入输出。因为其是 awk 的功能扩展，所以既可以独立以 bioawk 命令，也可以将其 link 到 awk/gawk ，直接代替 awk。

**内建函数**

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

**内建变量**

- bed: 1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts
- sam: 1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual
- vcf: 1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info
- gff: 1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute
- fastx: 1:name 2:seq 3:qual

### 1. 安装

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

### 2. 使用示例

使用示例部分使用 NCBI SRA 数据库的 illumina 测序数据 SRR1175124 为例

- fastq: SRR1175124
- bed:
- vcf:

#### 处理 fastq

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

#### 处理 bam

#### 处理 vcf

#### 处理 bed

#### 处理 gff
