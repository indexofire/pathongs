# Blast 使用说明

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    Blast 工具可以进行序列相似性比对，在 NGS 数据分析中经常会被使用到，特别是一些工具中需要 Blast 或 Blast+ 程序来作为第三方比对工具调用。拿到测序完成的草图后，因为基因组数据较大，连接NCBI网站往往又非常缓慢，所以要做 Blast 比对的话都需要做 Local Blast。这里介绍2个方式来实现：

    1. 用命令行进行本地 Blast
    2. 快速构建 Blast Web 服务

## 1. Blast 原理

## 2. 在线 blast

## 3. 基于命令行的 Local Blast

`Blast` 和 `NCBI Blast` 是 NCBI 先后发布的序列比对工具。`**B**asic **L**ocal **A**lignment **S**earch **T**ool`工具最早公布与1989年，我们一般都称呼其为`Blast`。2009年 NCBI 为了克服`Blast`的不足之处，为以后开发方便，重写代码，新的`NCBI Blast`（也称其为`Blast+`）工具应运而生。新版`Blast`工具在速度上有了提升，在输入输出上也更为灵活。

`[Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/)`工具最新版是2012年发布的**2.2.26**，目前已不再提供更新。`NCBI Blast`目前最新版为**2.9.0**。要区分2者也很简单，2个版本的程序运行方式不同：`Blast` 是通过 blastall -p 的方式调用子程序来比对搜索的，而`Blast+`则是直接使用blastn或blastp等子程序来比对搜索。另外前者用formatdb程序来格式化数据库，后者用makeblastdb程序来格式化数据。

### 3.1 Blast

**Blast的下载与安装**

```bash
# Ubuntu 中安装 blast
$ sudo apt install blast2

# 或者直接下载NCBI Linux预编译包，并解压缩安装
$ wget ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-x64-linux.tar.gz
$ tar -zxvf blast-2.2.26-x64 -C ~/app
# 在系统可执行路径中添加 blast 链接
$ sudo ln -s ~/app/blast-2.2.26-x64/blastall /usr/local/sbin

# 通过 conda 安装 blast
(blast)$ conda install blast-legacy
```

**运行 Blast**

`Blast`通过调用 blastall 程序来调用不同子程序实现不同类型序列比对。在命令行中输入`blastall`，会打印出一份参数列表。你也可以使用`man blast`来查看`Blast`工具的用户手册。

```bash
# 运行主命令
$ blastall
```

**参数说明：**

```bash
-p: 指定要运行的 blast 程序。可以调用的有:
    `blastp`             将氨基酸序列与蛋白质数据库比对
    `blastn`             将核酸序列与核酸数据库比对
    `blastx`             将核酸序列翻译成蛋白序列后与非冗余(nr)蛋白质数据库比对
    `psitblastn`
    `tblastn`
    `tblastx`
-d: 指定要调用数据库，默认值是非冗余(nr)数据库。本地 Blast 比对一般来说都是调用 formatdb 格式化的数据库。
-i: 输入，默认值是终端输入，也可以使用文件的方式，比如`-i seq.fasta`
-o: 输出，默认值是终端打印输出，也可以使用文件的方式，比如`-o result.txt`
-e: 期望值。
-T: 输出文件格式，默认值为F(False)，想输出HTML格式的可以用`-T T`
-M: 矩阵算法选择，默认值是`BLOSUM62`
-n: 如果想使用MegaBlast，就设置`-n T`
-b: 数据库中的比对结果显示条目数，默认是250条记录。
-m: 比对结果的输出显示方式，值为0~11，默认是0。各个数字含义见下表。
```

| 参数 | 含义 |
| ---- | ---- |
| 0 | pairwise |
| 1 | query-anchored showing identities |
| 2 | query-anchored no identities |
| 3 | flat query-anchored, show identities |
| 4 | flat query-anchored, no identities |
| 5 | query-anchored no identities and blunt ends |
| 6 | flat query-anchored, no identities and blunt ends |
| 7 | XML Blast output |
| 8 | tabular |
| 9 | tabular with comment lines |
| 10 | ASN, text |
| 11 | ASN, binary [Integer] |

其他参数参见`man blast`。

**应用举例:**

Local blast例子：

```bash
# 使用 custom_genome.fasta 文件作为数据库
# formatdb 格式化作为本地数据库
$ formatdb -i custom_genome.fasta -o T -p F
# 使用 blastn 对 myseq.fasta 序列进行比对
$ blastall -i myseq.fasta -d custom_genome.fasta -p blastn
```

### 3.2 Blast+

**下载安装 NCBI Blast+**

```bash
# Ubuntu 系统中安装 ncbi blast
$ sudo apt install ncbi-blast+

# Archlinux 系统中安装 ncbi blast
$ sudo pacman -S
```

**构建数据库**

```bash
# 将序列文件 data/database.fasta 格式化成数据库
# 数据类型为核酸，则参数为 -dbtype nucl
$ makeblastdb -in data/database.fasta -dbtype nucl -parse_seqids
```

**运行 blast**

```bash
# 使用 blastn 程序比对 myseq.fasta 序列。
$ blastn -query myseq.fasta -db data/database.fasta -out result.txt
# 常用简单输出方式
$ blastn -query myseq.fasta -db data/database.fasta -outfmt 6 -out result.txt
```

http://www.personal.psu.edu/iua1/courses/files/2014/lecture-12.pdf

## 4. 构建本地的 blast web 服务

### 4.1 blastkit

blastkit 是一个包含webserver等工具的blast工具集。

**安装依赖包**

```bash
$ conda create -n blastkit
$ conda activate blastkit
(blastkit)$ pip install pygr
(blastkit)$ pip install whoosh
(blastkit)$ pip install git+https://github.com/ctb/pygr-draw.git
(blastkit)$ pip install git+https://github.com/ged-lab/screed.git
(blastkit)$ conda install lighttpd
(blastkit)$ conda install blast-legacy
```

**对lighttpd webserver进行配置**

```bash
$ cd /etc/lighttpd/conf-enabled
$ sudo ln -fs ../conf-available/10-cgi.conf ./
$ sudo echo 'cgi.assign = ( ".cgi" => "" )' >> 10-cgi.conf
$ sudo echo 'index-file.names += ( "index.cgi" ) ' >> 10-cgi.conf
$ sudo /etc/init.d/lighttpd restart
```

```bash
(prokka)$ conda install lighttpd
(prokka)$ cd $HOME/.conda/etc/lighttpd/conf-enabled
$ sudo ln -fs ../conf-available/10-cgi.conf ./
$ sudo echo 'cgi.assign = ( ".cgi" => "" )' >> 10-cgi.conf
$ sudo echo 'index-file.names += ( "index.cgi" ) ' >> 10-cgi.conf
$ lighttpd restart
```


**安装 blastkit**

```bash
# 克隆软件仓库 blastkit
$ git clone https://github.com/ctb/blastkit.git -b ec2
$ cd blastkit/www
# 将www路径链接到 /var 中
$ sudo ln -fs $PWD /var/www/blastkit
# 建立临时文件存放目录
$ mkdir files
$ chmod a+rxwt files
#
$ chmod +x /root
$ cd ..
# 运行检测脚本
$ python ./check.py

#如果安装顺利，就会提示一切已经准备完毕。接下来要准备数据。
# 将大肠埃希菌DNA序列文件contigs.fa复制到数据库目录下
$ cp ~/contigs.fa blastkit/db/db.fa
# 使用 formatdb 格式化数据库
$ formatdb -i db/db.fa -o T -p F
# 建立索引
$ python index-db.py db/db.fa
```

## Reference

* [Blastkit](https://github.com/ctb/blastkit.git)
* [Caltech workshop](https://github.com/dib-lab/2013-caltech-workshop/blob/master/blastkit.txt)
* [Angus](http://angus.readthedocs.io/en/2017/running-command-line-blast.html)
