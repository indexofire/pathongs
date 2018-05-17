# Local Blast 使用说明

> Blast 工具可以进行序列相似性比对，在 NGS 数据分析中经常会被使用到，特别是一些工具中需要 Blast 或 Blast+ 程序来作为第三方比对工具调用。拿到测序完成的草图后，因为基因组数据较大，连接NCBI网站往往又非常缓慢，所以要做 Blast 比对的话都需要做 Local Blast。这里介绍2个方式来实现：

1. 用命令行进行本地 Blast
2. 快速构建 Blast Web 服务

---

## 1. 基于命令行的 Local Blast

`blast` 和 `blast+` 这2个程序集容易搞混。NCBI 最早在1989年创建`Basic Local Alignment Search Tool`工具，沿用至2009年无论是命令行工具或是在线程序，都称呼其为`blast`。2009年 NCBI 鉴于 `blast` 的一些不足，重新开发了新的`blast+`命令行工具，新的 blast+ 工具在速度上有了提升，在输入输出上也更为灵活。

目前 Blast 工具最新版是2012年发布的**2.2.26**（可能已经处于暂停升级的状态？），而 Blast+ 目前一直在更新。要区分2者也很简单，blast 是通过 blastall -p 的方式调用子程序来比对搜索的，而 blast+ 则是直接使用 blastn 或 blastp 等子程序来比对搜索。另外前者用 formatdb 程序来格式化数据库，后者用 makeblastdb 程序来格式化数据。

### 1.1 Blast

**Blast的下载与安装**

```bash
# 安装 Ubuntu 编译包。
$ sudo apt-get install blast2

# 或者直接下载NCBI Linux预编译包，并解压缩安装
$ wget ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-x64-linux.tar.gz
$ tar -zxvf blast-2.2.26-x64 -C ~/app
$ sudo ln -s ~/app/blast-2.2.26-x64/blastall /usr/local/sbin
```

**运行 Blast**

Blast 通过调用 blastall 这个 gateway 程序，来分别调用不同算法和程序实现序列比对。在命令行中输入`blastall`，会打印出一份参数列表。你也可以使用 `man blast` 来查看 Blast 工具的用户手册。

```bash
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

| options | meaning |
| -- | -- |
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

Local blast例子：首先下载一个基因组文件并格式化作为本地数据库，然后使用 blastn 对序列进行比对。

```bash
$ formatdb -i custom_genome.fasta -o T -p F
$ blastall -i myseq.fasta -d custom_genome.fasta -p blastn
```

### 1.2 Blast+

**下载安装 NCBI Blast+**

```bash
$ sudo apt-get install ncbi-blast+
```

**构建数据库**

```bash
$ makeblastdb -in data/database.fasta -dbtype nucl -parse_seqids
```

**运行 blast**

```bash
$ blastn -query my_seq.fasta -db data/database.fasta -out result.txt
```

http://www.personal.psu.edu/iua1/courses/files/2014/lecture-12.pdf

## 2. 构建本地的 blast web 服务

### 2.1 blastkit

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
$ cd ~/app
$ git clone https://github.com/ctb/blastkit.git -b ec2
$ cd blastkit/www
~/app/blastkit/www$ sudo ln -fs $PWD /var/www/blastkit
~/app/blastkit/www$ mkdir files
~/app/blastkit/www$ chmod a+rxwt files
~/app/blastkit/www$ chmod +x /root
~/app/blastkit/www$ cd ..
~/app/blastkit$ python ./check.py
```

如果安装顺利，就会提示一切已经准备完毕。接下来要准备数据。

```bash
$ cp /mnt/assembly/ecoli.21/contigs.fa ~/app/blastkit/db/db.fa
$ cd ~/app/blastkit
$ formatdb -i db/db.fa -o T -p F
$ python index-db.py db/db.fa
```

## Reference

* [Blastkit](https://github.com/ctb/blastkit.git)
* [Caltech workshop](https://github.com/dib-lab/2013-caltech-workshop/blob/master/blastkit.txt)
* [Angus](http://angus.readthedocs.io/en/2017/running-command-line-blast.html)
