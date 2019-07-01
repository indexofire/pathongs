# Prokka

---

[Prokka][] 是一个原核物种基因组注释工具，由墨尔本大学生物信息学家 [Torsten Seemann](https://tseemann.github.io/) 开发的基于命令行的本地快速注释工具，主要用来注释小基因组比如细菌、病毒等。

## 安装

### 1. 编译包安装

如果下载预编译包安装的话，需要有Bioperl支持。以及一些第三方库和软件的支持。

**安装 Bioperl**

```bash
$ sudo apt-get install cpanminus
$ cpanm Bio::Perl XML::Simple
```

**安装其他依赖**

```bash
$ sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
```

**安装 rRNA 注释功能**

这里选择 Barrnap 和 MINCED

```bash
$ wget http://www.vicbioinformatics.com/barrnap-0.6.tar.gz
$ tar zxvf barrnap-0.5.tar.gz -C ~/app/
$ echo "export PATH=$PATH:$HOME/app/barrnap-0.6/bin" >> ~/.bashrc
$ source ~/.bashrc

$ wget https://github.com/ctSkennerton/minced/releases/download/0.3.0/minced
$ sudo cp minced /usr/local/sbin
```

**安装 ncRNA 注释功能**

使用 infernal

```bash
$ wget http://eddylab.org/infernal/infernal-1.1.2.tar.gz
$ tar zxvf infernal-1.1.2.tar.gz -C ~/app/
$ echo "export PATH=$PATH:$HOME/app/infernal-1.1.2"
```

**安装 prokka**

```bash
$ wget http://www.vicbioinformatics.com/prokka-1.12.tar.gz
$ tar zxvf prokka-1.12.tar.gz -C ~/app
$ echo "export PATH=$PATH:$HOME/app/prokka-1.12/bin" >> ~/.bashrc
$ source ~/.bashrc
$ prokka --setupdb
```

如果要使用 --gram 功能，需要安装 [signalp](http://www.cbs.dtu.dk/services/SignalP/)。[signalp](http://www.cbs.dtu.dk/services/SignalP/) 是 [DTU 大学生物及卫生信息学系](http://www.cbs.dtu.dk/index.html)开发，支持web界面的工具。命令行版本需要提供非商业ISP提供的邮箱帐号注册以后才能下载（比如政府或大学）。点击[这里](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)。

### 2. conda 安装

```bash
$ conda create -n prokka
$ conda activate prokka
(prokka)$ conda install prokka
```

## 使用 Prokka

要想注释fasta文件，最简单的方式，不填加任何参数就可以获得。所有的结果会放在

```bash
# 下载基因组数据
$ esearch -db nuccore -query "1639[txid] AND egd" | efetch -format fasta > egd.fasta

# 注释数据
$ prokka egd.fasta
```

通过添加`--listdb`可以查看数据库信息。

```bash
# 查看可以使用的数据库
$ prokka --listdb

Looking for databases in: ...
* Kingdoms: Archaea Bacteria Mitochondria Viruses
* Genera: Enterococcus Escherichia Staphylococcus
* HMMs: HAMAP
* CMs: Bacteria Viruses
```

如果需要比较完整的进行序列注释，结果可以比较方便的用来提交的话，可以添加相应的参数来规范注释以及调用不同的工具扫描ncRNA等。

```bash
# 注释 Listeria monocytogenes 标准株 egd
$ prokka --outdir egd --prefix egd --addgenes --addmrna --compliant --centre CDC --genus Listeria --species "Listeria monocytogenes" --strain egd
--kingdom Bacteria --usegenus --cpus 4 --rfam --rnammer --force egd.fasta

# 生成结果注释文件
$ ls -l egd | awk '{print $9}'

egd.err     # prokka 对注释结果存在的一些疑问进行报告的信息
egd.faa     # 注释的氨基酸序列
egd.ffn     # 注释的碱基序列
egd.fna     # 以 NCBI gnl|centre| 为ID命名的碱基序列文件
egd.fsa     # 以 NCBI gnl|centre| 为ID命名的碱基序列文件
egd.gbk     # genbank 格式的注释文件
egd.gff     # gff 格式的注释文件
egd.log     # prokka 运行日志
egd.sqn     # sqn 格式的文件，可以用来提交到 NCBI
egd.tbl     # tbl 格式的文件，可以用来提交到 NCBI
egd.tsv     # tsv 格式文件，
egd.txt     # prokka 注释的各种类型序列统计信息
```

**--kingdom**

默认注释的是细菌基因组，如果是其他物种则要添加这个参数。可选项有：

* Archaea
* Bacteria
* Mitochondria
* Viruses

比如要注释病毒基因组：

```bash
$ prokka --kingdom Viruses contigs.fasta
```

### .err 文件


### .tsv 文件

tsv 文件按照 locus_tag 顺序排序了注释的结果。字段包括`locus_tag`, `ftype`, `length_bp`, `gene`, `EC_number`, `COG`, `product`。

```bash
$ head egd.tsv
locus_tag       ftype   length_bp       gene    EC_number       COG     product
DJECODEN_00001  CDS     1356    dnaA            COG0593 Chromosomal replication initiator protein DnaA
DJECODEN_00001  gene    1356    dnaA
DJECODEN_00001  mRNA    1356    dnaA
DJECODEN_00002  CDS     1146    dnaN            COG0592 Beta sliding clamp
DJECODEN_00002  gene    1146    dnaN
DJECODEN_00002  mRNA    1146    dnaN
DJECODEN_00003  CDS     1344    yeeO_1          COG0534 putative FMN/FAD exporter YeeO
DJECODEN_00003  gene    1344    yeeO_1
DJECODEN_00003  mRNA    1344    yeeO_1
```

```bash
# 注释CDS长度排序
$ awk '/CDS/ {print $0}' egd.tsv | sort -n -k3 -r | less

# 显示所有的不重复注释产物
$ awk '!a[$1]++{print}' egd.tsv

# 显示所有COG注释
$ awk '/COG/' egd.tsv
```

[Prokka]: https://github.com/tseemann/prokka "Prokka"
