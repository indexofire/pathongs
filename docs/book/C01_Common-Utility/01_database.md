# 基因组公共数据库

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    测序的数据如原始 reads 数据我们可以上传到 NCBI 或 EBI 的公共数据库，发表的文章中可以使用。此外公共数据库中的数据也可以为我们的研究提供帮助。本节我们就来了解一下如何上传和下载基因组测序数据。首先让我们来了解一下网络上的测序公共数据库：

**常见的高通量测序数据的公共数据库：**

- [SRA][] 短序列数据库：由 [NCBI][] 负责维护
- [ENA][] 欧洲核酸数据库：由 [EBI][] 负责维护
- [GSA][] 中国组学数据库：由[中科院北京基因研究所](http://www.big.ac.cn)负责维护

[SRA][] 是 [NCBI][] 为了应对越来越多的高通量测序数据而在 2007 年底推出的测序数据库，用于存储、显示、提取和分析高通量测序数据。而 [ENA][] 则是由 [EBI][] 负责维护的功能类似的数据库，同时作为 [Ensembl](http://www.ensembl.org)、[UniProt](http://www.uniprot.org) 和 [ArrayExpress](http://www.ebi.ac.uk/arrayexpress) 等服务的底层基础。2者在主要功能方面非常类似，同时数据互通。

---

## 1.SRA 数据库

### 1.1 简介

[SRA][] 是 Sequence Read Archive 的首字母缩写。[SRA][] 与 Trace 最大的区别是将实验数据与 metadata（元数据）分离。metadata 是指与测序实验及其实验样品相关的数据，如实验目的、实验设计、测序平台、样本数据(物种，菌株，个体表型等)。metadata可以分为以下几类：

* Study：accession number 以 DRP，SRP，ERP 开头，表示的是一个特定目的的研究课题，可以包含多个研究机构和研究类型等。study 包含了项目的所有 metadata，并有一个 NCBI 和 EBI 共同承认的项目编号（universal project id），一个 study 可以包含多个实验（experiment）。
* Sample：accession number以 DRS，SRS，ERS 开头，表示的是样品信息。样本信息可以包括物种信息、菌株(品系) 信息、家系信息、表型数据、临床数据,组织类型等。可以通过[Trace](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_obj) 来查询。
* Experiment：accession number 以 DRX，SRX，ERX 开头。表示一个实验记载的实验设计（Design），实验平台（Platform）和结果处理（processing）三部分信息。实验是 [SRA][] 数据库的最基本单元，一个实验信息可以同时包含多个结果集（run）。
* Run：accession number 以DRR，SRR，ERR 开头。一个 Run 包括测序序列及质量数据。
* Submission：一个 study 的数据，可以分多次递交至 SRA 数据库。比如在一个项目启动前期，就可以把 study，experiment 的数据递交上去，随着项目的进展，逐批递交 run 数据。study 等同于项目，submission 等同于批次的概念。

### 1.2 数据下载

这些公共数据库的数据可以通过不同方式下载，如https、ftp，或者使用aspera这个工具进行高速下载。

#### 1.2.1 用 Aspera 工具下载 SRA

**Aspera 介绍**

[Aspera][] 是 IBM 开发的加速在线数据传输的一组服务器端与客户端的应用程序。[Aspera][] 的各个工具可以在 [官方网站](http://downloads.asperasoft.com/downloads) 下载获得。

在服务器上我们只需要命令行工具 `ascp` 即可。而对于桌面用户来说，可以选择安装 [Aspera Connent](http://downloads.asperasoft.com/en/downloads/8?list) 等工具。由于 Ascp Client 工具需要登录后才能下载，因此选择了 [Aspera Connect][] 来安装。它是 [Aspera][] 用于浏览器下载的高效插件，其内嵌了 `ascp` 命令也可以直接在服务器上运行。安装了 [Aspera Connect][] 的带图形界面的客户端，浏览器会在有 Connect, Faspex 或者 Shares 页面的链接处添加 ascp 下载快捷按钮。

**安装 aspera connect**

[Aspera Connect][] 下载解压缩后，可以看到一个 .sh 文件，执行后运行后会安装到当前用户 $HOME.aspera/connect 目录下。

```bash
# 下载并解压缩程序
$ wget http://download.asperasoft.com/download/sw/connect/3.6.1/aspera-connect-3.6.1.110647-linux-64.tar.gz
$ tar zxf asper-connect-3.6.1.110647-linux.tar.gz

# 安装到本地用户
$ ./aspera-connect-3.6.1.110647-linux-64.sh

# 添加用户执行路径
$ echo 'PATH="$PATH:$HOME/.aspera/connect/bin/"' >> ~/.bashrc
```

对于喜欢使用客户端的桌面用户，也可以使用 [Aspera Desktop Client][] 客户端程序来下载数据。

```bash
$ sed 's/^INSTALL_DIR/INSTALL_DIR=~/envs/srst2/src/ascp/connect/' aspera-connect-3.6.1.110647-linux-64.sh > install.sh
$ ./install.sh
```

**下载 SRA 数据**

!!! warning "注意"
    注意新版的`ascp`用`.openssh`作为密钥文件而不是原来的`.putty`

```bash
# ascp 下载 SRA 数据
$ ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
> --user=anonftp --host=ftp.ncbi.nlm.nih.gov --mode=recv -l 100m -pQTk1 \
> /sra/sra-instant/reads/ByRun/sra/SRR/SRR955/SRR955386/SRR955386.sra .

# 添加 alias 简化操作
$ echo "alias ascp='ascp -i ~/.ascpera/connect/etc/asperaweb_id_dsa.openssh --user=dbtest --host=sra-download.ncbi.nlm.nih.gov --mode=recv -l 100m -pQTk1'" >> ~/.bashrc
$ source ~/.bashrc

# 测试下载
$ ascp /data/sracloud/srapub/ERR579922 ./ERR578822.sra
```

#### 1.2.2 ftp 下载

从国内公网访问 NCBI FTP 服务器往往是非常缓慢的，而 [EBI][] 就不同了，有些地方的速度很快，下载甚至可以超过 [Aspera][]。由于 [NCBI][] 和 [EBI][] 共享数据库，因此对于国内 [NCBI][] 下载缓慢或者无法用 ascp 的用户，可以考虑在 [EBI][] 下载数据。

### 2.2 sratoolkit

[NCBI][] 开发了 [sratoolkit][] 工具来帮助处理 SRA 数据，正确配置后可以很方便的下载 SRA 数据。

#### 2.2.1 安装sratoolkit

可以直接从 [NCBI][] 上下载。最新的源码可以在 [Github](https://github.com/ncbi/sra-tools) 获得。

```bash
# download from NCBI
$ cd /tmp
$ wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.6/sratoolkit.2.5.6-ubuntu64.tar.gz
$ tar zxf sratoolkit.2.5.6-ubuntu64.tar.gz -C ~/app
```

#### 2.2.2 下载 SRA 格式文件

如果你安装并设置了 [Aspera Connect][]，那么`prefetch`会优先使用`ascp`方式来下载，如果没有安装或则`ascp`下载失败，则切换成 HTTP 方式下载 sra 数据。另外`fastq-dump`命令也能从远端直接下载数据，加上`-X 1`参数，会预先下载最前的5个 reads，加上`-Z`参数，则会将这些 reads 打印到终端输出。

```bash
# 下载 SRR955386.sra 到 你安装 sratoolkit 时配置的 public 目录中，默认是 ~/ncbi/public/sra
$ prefetch -v SRR955386
# 下载 SRR955386.sra 到 你安装 sratoolkit 时配置的 public 目录中，默认是 ~/ncbi/public/sra，并且在终端输出5行 reads 数据。
$ fastq-dump -X 5 -Z SRR955386
```

#### 2.2.3 并转换成 .fastq 格式

获得了 .sra 文件后，需要将其转换成 .fastq 格式的文件，用`fastq-dump`可以很方便的实现。转换之前要注意的是该 run 的 metadata 里，测序类型是 SE 还是 PE 的。

```bash
# 将 sra 文件移动到 ~/data 目录中
$ mv ~/.ncbi/public/sra/SRR955386.sra ~/data
# 如果是 SE 测序数据，直接运行以下命令
$ fastq-dump SRR955386.sra
# 如果是 PE 测序数据，则要添加参数：--split-files
$ fastq-dump --split-files SRR955386.sra
```

#### 2.2.4 SE/PE 文件判断

正常的 sra 文件的 metadata 应该包含测序采用的是 SE 还是 PE 的方式。但如果你不知道所下载的到底是 SE 还是 PE 格式的文件可以用`fastq-dump -X 1 --split-spot`的方法来判断。

```bash
# it's SE if nreads=1, and PE when nreads=2，统计整个文件，因此速度比较慢
$ sra-stat -xs SRR955386.sra | grep "nreads"

# 如果输出是4，那么就是SE，如果是8,则是PE
$ fastq-dump -X 1 --split-spot -Z SRR955386.sra | wc -l

# 或者加上参数让 fastq-dump 自己判断
# 当 sra 文件是 SE 测序时，fastq-dump只能dump出1个 *_1.fastq 文件
$ fastq-dump --split-files ERR493452.sra
$ mv ERR493452_1.fastq ERR493452.fastq
```

当需要判断批量下载的 sra 文件时，区分那些是 PE 的那些是 SE 的文件，可以用以下脚本：

``` python
import os
import subprocess


def check_SRA_type(fn):
    fn = os.path.abspath(fn);
    try:
        contents = subprocess.check_output(["fastq-dump", "-X", "1", "-Z", "--split-spot", fn]);
    except subprocess.CalledProcessError, e:
        raise Exception("Error running fastq-dump on", fn);
    # -X 1 will output 4 lines if SE, and 8 lines if PE
    if(contents.count("\n") == 4):
        return False;
    elif(contents.count("\n") == 8):
        return True:
    else:
        raise Exception("Unexpected output from fastq-dump on ", filename);
```

#### 2.2.5 利用entrez批量下载

如果想下载一个完整的 project 数据，可以利用 entrezdirect 工具。

```bash
~/tmp$ wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
~/tmp$ tar zxvf edirect.tar.gz -C ~/app
~/tmp$ sudo ln -s ~/app/edirect/* /usr/local/sbin/
~/data$ esearch -db sra -query PRJNA40075  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files
```

如果想下载不同 research 的数据，可以自己建立一个 accession number list 的文件，比如利用 [NCBI] 的 entrez 网页界面，导出所需要的数据 accession number，然后利用 `ascp` 下载。不过建议`ascp`下载不要太狠心，否则容易被 [NCBI][] 给封掉。

## 3. Assembly 数据库

如果不需要研究 reads 数据，而是只想研究拼接的基因组数据，可以从 NCBI assembly 数据库去寻找数据，并通过 ascp/ftp 等工具下载。

## 4. 上传数据

SRA 上传测序结果可以参照 [NCBI文档](http://www.ncbi.nlm.nih.gov/books/NBK47529/) 来实现。一般上传数据到NCBI SRA的过程需要6步：

1. Create a BioProject for this research
2. Create a BioSample submission for your biological sample(s)
3. Gather Sequence Data Files
4. Enter Metadata on SRA website
  4.1 Create SRA submission
  4.2 Create Experiment(s) and link to BioProject and BioSample
  4.3 Create Run(s)
5. Transfer Data files to SRA
6. Update Submission with PubMed links, Release Date, or Metadata Changes

需要注意的一点是，上传的过程中很多地方一旦保存或提交就不可以修改，尤其是各处的Alias。但是，可以联系NCBI的工作人员修改内容。NCBI的工作效率是很高的，一般不超过48小时，就可以得到确认，并拿到登录号。


## [BIG(中国科学院北京基因组所)](http://bigd.big.ac.cn)

几个数据库入口：

- [Bioproject](http://bigd.big.ac.cn/bioproject)
- [Biosample](http://bigd.big.ac.cn/biosample)
- [GSA](http://bigd.big.ac.cn/gsa)

一般可以通过FTP方式下载:w


## [复制起点数据库](http://tubic.org/doric)

细菌、真菌和质粒的复制起始位点数据库


## [质粒数据库](http://www.patlas.site)

## [plsdb质粒数据库](http://ccb-microbes.cs.unisaarland.de/plsdb)

## [eggNOG](http://eggnog.embl.de)

## [GO](http://geneontology.org)


## Reference

1. [NCBI上传数据文档](http://www.ncbi.nlm.nih.gov/books/NBK47529/)
2. 熊筱晶, NCBI高通量测序数据库SRA介绍, 生命的化学[J], 2010:6, 959-963.
3. http://blog.sciencenet.cn/blog-656335-908140.html
4. https://www.biostars.org/p/139422/
5. https://www.youtube.com/watch?v=NSIkUHKRPpo

[NCBI]: http://www.ncbi.nlm.nih.gov
[SRA]: http://www.ncbi.nlm.nih.gov/sra
[EBI]: http://www.ebi.ac.uk
[ENA]: http://www.ebi.ac.uk/ena
[Aspera]: http://asperasoft.com
[Aspera Connect]: http://download.asperasoft.com/download/docs/connect/3.6.1/user_linux/webhelp/index.html#dita/introduction.html
[Aspera Desktop Client]: http://downloads.asperasoft.com/en/downloads/2
[GSA]: http://gsa.big.ac.cn/
[sratoolkit]: (https://github.com/ncbi/sra-tools)
