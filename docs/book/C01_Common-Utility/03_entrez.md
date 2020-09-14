# Entrez

{{ git_page_authors }} 更新于: {{ git_revision_date }}

![edirect](../../assets/images/C02/03/banner.jpg)

---

## Entrez 是什么

---

## 玩转 entrez-direct

!!! abstract "内容简介"
    [NCBI][] 的 [Entrez][] 工具功能非常强大，既可以从 web 页面访问 [NCBI][] 来查询，也可以利用 [Entrez][] 提供的 [Web Services](https://www.ncbi.nlm.nih.gov/books/NBK25501/) 来实现诸多功能（编写第三方程序等）。此外 [NCBI][] 还提供了 [Entrez][] 的命令行工具 edirect[^1]。edirect 全名为 **Entrez Direct**，里面包含了一组各司其职的工具和脚本，通过 Linux 系统下的管道pipeline[^4]功能联合使用这些命令行工具，在服务器上实现高效的进行 [Entrez][] 的检索，抓取，过滤，排序等操作。

### 1. 安装

#### 1.1 下载安装

edirect 直接下载预编译包，添加到系统路径中即可。

```bash
# 下载 edirect 安装包
$ wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz

# 将 edirect 安装在用户 app 目录下，可以根据需求自行选择安装路径
$ tar zxvf edirect.tar.gz -C ~/apps

# 运行安装脚本
$ cd ~/apps/edirect
$ bash ./setup.sh

# 添加软件安装目录到系统环境变量中
$ echo "export PATH='\$PATH:\$HOME/apps/edirect/'" >> ~/.bashrc
$ source ~/.bashrc
```

!!! note "安装说明"
    耐心等待一会，当看到下面字符时，表示安装成功。
    ```
    ENTREZ DIRECT HAS BEEN SUCCESSFULLY INSTALLED AND CONFIGURED
    ```

也可以不运行 `setup.sh` 直接使用 edirect 的一些工具，不过如果要使用全部功能，有些perl模块需要安装，运行 `setup.sh` 来自动完成 CPAN 的 modual 安装以及 go 版本的 xtract 的编译（不运行安装脚本会调用 perl 版本的 xtract，数据量大时速度略慢）。

!!! hint "提示"
    对于部分操作系统如 **ubuntu** 等，**uname** 位于 /bin，而 xtract 接口调用脚本使用 /usr/bin/uname 来判断使用那个版本的 xtract，如果要在这些操作系统下正常使用 xtract，需要让程序在正确路径可识别，这里通过添加系统 link 到 /usr/bin/unname。

```bash
# 查看 uname 路径。
$ which uname

# 如果是位于 /bin/ 下，将其链接到 /usr/bin/ 下即可。
$ sudo ln -s /bin/uname /usr/bin/uname
```

#### 1.2 conda 中安装

如果系统通过 conda 进行软件管理，可以通过如下方式安装。

```bash
$ conda install entrez-direct
```

不过 conda 安装的还是 perl 版本的 xtract，而不是 go 版本的，需要到 conda/pkgs 包根目录中的 edirect 相应目录，用 `setup.sh` 自行编译后生成 xtract.Linux(Linux 操作系统)，将其复制到 conda 的 bin 路径即可。go 版本的 xtract 可以用 -help 参数查看帮助。

### 2. 基本用法

表. 最常用的工具参见下表

| 工具名称 | 工具用途 | 常用参数 |
| -------- | -------- | ---- |
| **esearch** | 搜索命令，将所要检索的内容提交到 Entrez 中，返回相应的结果记录 | -db、-query |
| **efetch** | 下载 NCBI 数据库中的记录和报告并以相应格式打印输出 | -db、-id、-format、-mode |
| **einfo** | 获取目标结果在数据库中的信息 | -db、-dbs、-fields、-links |
| **elink** | 对目标结果在其他数据库中比配结果 | -db、-id、-related、-target、-name |
| **epost** | 上传 UIDs 或者 序列登记号 | -db、-id、-format、-input、-label|
| **efilter** | 对之前的检索结果进行过滤或限制 | -query、-sort、-field |
| **xtract** | 将esearch获得的 XML 格式结果转换成表格格式 | -pattern、-if、-block、-element、-sep、-filter |
| **esummary** | 获得 XML 格式的建立 | -db、-id、-format、-mode |
| **ecitmatch** | 统计引用数据 |  -journal、-year、-volume、-page、-author |

#### 2.1 esearch

**esearch** 负责检索并返回 xml 格式的命中情况。使用 esearch 命令，相当于平时你通过浏览器访问 NCBI，在 [GQuery](https://www.ncbi.nlm.nih.gov/gquery/) 的 [nucletide](https://www.ncbi.nlm.nih.gov/nuccore/) 数据库搜索见面搜素关键字 `salmonella`，返回命中结果。

```bash
$ esearch -db nuccore -query "salmonella"
```

终端里会返回命中结果：

``` xml
<ENTREZ_DIRECT>
  <Db>nuccore</Db>
  <WebEnv>NCID_1_153624467_130.14.18.34_9001_1490828308_1411373464_0MetA0_S_MegaStore_F_1</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>2093561</Count>
  <Step>1</Step>
</ENTREZ_DIRECT>
```

**esearch** 需要设置的参数有 **-db**（所要搜索的目标数据库）和 **-query**（所要搜索的关键词）。

结果默认返回 xml 格式的数据（参见 ```<ENTREZ_DIRECT>...</ENTREZ_DIRECT>```），其中 **<Count>** 中的数值表示检索命中的记录数，这与 [NCBI][] 的 [Entrez][] 网页界面里进行检索操作的结果是一致的。

**esearch** 还有一个参数 -sort 可以根据数据库类型返回不同的排序方式。默认的排序方式是按照时间逆序排列的。想要得到准确的查询结果，优化 query 很重，可以参见 [Filter](#3-filtering-queries) 相关内容。

#### 2.2 efetch

**efetch** 可以用来下载各个数据库的各种格式数据。你可以通过2种方式下载：

1. 用 esearch 获得命中记录后再下载。下面这个命令，就会下载所有结果并保存到 result.fasta 文件中。

```bash
$ esearch -db nuccore -query "salmonella" -days 10 | \
> efetch -format fasta -db nuccore > result.fasta
```

!!! warning
    一定要注意搜索命中值的数量，关键词命中太多，数据量会太大，全部下载可能直接进入假死状态。应该通过设置参数来限定搜索结果。当不确定时，先用esearch工具看返回的命中值<Count></Count>有多少再决定是否进一步限定。

2. 直接指定 id 下载

```bash
$ efetch -db pubmed -id 2142131 -format txt
```

**efetch** 设置的参数是 **-db**，将搜索命中的记录从该数据库进行。如果 **esearch** 检索的结果和需要下载的结果来自同一数据库，可以不添加 **-db** 参数，系统会默认使用esearch的数据库。

**efetch** 还要设置的参数是 **-format**，即下载数据的格式。常用的有 **xml** 格式，后续可以通过 **xtract** 和 **linux** 一些命令行工具来做处理。你也可以直接添加 **-id** 参数直接进行下载而不通过 **esearch** 检索后下载。如下面例子中将 **pubmed id 2137412** 的文献输出，格式为 Medline[^2]

```bash
# 搜索 pubmed 数据库编号 2137412 的文献，以 medline 格式输出
$ efetch -db pubmed -id 2137412 -format medline
```

对于 efetch 可以下载那些数据库，对于各个数据库名称和数据格式，可以用`-h`来查看

```bash
$ efetch -h
```

#### 2.3 xtract

edirect 的核心命令，由于 **esearch** 查询和下载的数据常常是 xml 格式的，对于 XML[^3] 格式的数据，可以用 **xtract** 工具将其转换成纯文本格式的内容。xtract 不仅可以转换 [NCBI][] 的 XML 格式文件，对于任何通用的 XML 格式文件都可以进行转换。

xtract 基本用法：

```
# 读取 XML 文件方式
$ xtract -input XML-file.xml
```

一般我们都是通过管道工具 Pipeline 组合 esearch 和 efetch 等工具获得 XML 格式的 NCBI 数据库命中记录后，再用 xtract 进行处理。

```bash
$ esearch -db pubmed -query 'E. coli AND (2016/12/30[pdat]:2016/12/31[pdat])' | \
> efetch -format xml | xtract -pattern PubmedArticle -element \
> MedlineCitation/PMID Journal/ISOAbbreviation ArticleTitle
```

这条命令将检索到的文献记录转换成行，再以 PMID，杂志名缩写和文章题目为列展开。命令打印的结果如下：

```
28239441J Public Health AfrQuality and labeling information of Mo       ringa oleifera products marketed for HIV-infected people in Zimbabwe.
28115890Korean J Food Sci Anim ResourAntimicrobial Activity of Kefir against Variou     s Food Pathogens and Spoilage Bacteria.
28110668Enzyme Microb. Technol  .A novel thermophilic and halophilic esterase from Janibacter sp. R02,   the first member of a new lipase family (Family XVII).
...
27030192J.   Vet. Sci.Development and evaluation of an immunochromatographic assay   using a gp51 monoclonal antibody for the detection of antibodies against the bovine leukemia virus.
```

当我们不清楚XML元素结构时，可以用以下命令打印出树形结构。

```bash
$ xtract -input my.xml -outline
```

xtract 可以灵活的改变数据的格式，所以涉及的参数也比较多。要想获得自己需要的数据格式，要想用好 edirect，掌握好 xtract 的用法是必不可少的。

**-pattern 参数**: 将 xml 格式数据，将<xml element>转换成行格式数据。每一个element包含的数据做为一行。

xtract 命令首先要通过 -pattern 定义所要分析的数据级别。比如对于 pubmed 数据库，efetch 输出 -format 格式为 docsum，那么我们一般要抓取的是每一个命中的文献记录，-pattern 就可以用 DocumentSummary 元素。如果你只想知道所有命中文献的作者（比如按照次数排序），那么就可以直接将 -pattern 定义为 Author。

**-element 参数**: 将 **-pattern** 定义的行中以每个element形成一列。每个element以tab键分隔。

每个 -pattern 定义后，要多 -pattern 定义的元素的子元素进行输出，就要用到 -element 参数了。-element 可以跟多个子元素，默认以tab键分割。但是要注意的是子元素必须是唯一子元素，比如 -pattern 如果定义为 DocumentSummary，那么每个 DocumentSummary 元素有一个 Authors 子元素，而每个 Authors 有多个 Author 子元素（因为一篇文献有多个作者），如果直接将 -element Author，是无法返回作者信息的。需要用 -block 参数把 Authors 列转换成行再进行分割。

**-block 参数**: 可以将某个XML标签下的多个相同标签转换成一行的多个元素，最常见的用法是在 `pubmed` 数据库检索文献时<AuthorList>下多个<Author>标签打印出来。

```bash
$ esearch -db pubmed -query "Salmonella[ORGN]" -days 5 | efetch -format xml | \
> xtract -pattern PubmedArticle -block Author -sep " " -element LastName,Initials
```

**-tab、-sep 参数**: 将多条 element 记录分割。

**/**: 主要用在定义连续多层标签的定义上。比如PMID标签即可以在 MedlineCitation 标签下，有时也会在 RefSource 标签下。所以要明确检索文献的 PMID，就应该用 **MedlineCitation/PMID**

**@**: 指定 XML 元素 element 的 attribute 值

**-if、-if/-equals、-if/-contains 参数**: 条件运算符

**-and、-or、-position、-def 参数**: 关系运算符

#### 2.4 elink

当我们要检索一个数据库数据所关联到其他数据库时，就要用到elink工具了。比如我们先在SRA数据库中搜索到所有鼠伤寒沙门菌的SRA实验数据，然后链接到biosample数据库，获得这些菌株的生物样本信息。

```bash
$ esearch -db sra -query "salmonella typhimurium" | elink -target biosample | efetch -format docsum > result.xml
```

#### 2.4 einfo

对于 NCBI 数据库的各项复杂参数，字段等要完全了解和掌握比较困难。edirect 里提供了 einfo 工具，这样其他工具所能访问的具体数据库、字段和可链接数据库等信息都可以用 einfo 了解。

```bash
# 查询 `esearch -db` 参数所能访问的所有数据库名称
$ einfo -dbs | sort

# 查询某个数据库字段情况，如 sra 数据库所有对应的字段名称
$ einfo -db sra -fields

# 查询数据库所有可链接的其他数据库
$ einfo -db sra -links
```

### 3. 高级检索

工具 esearch 和 efilter 在做检索时如果做好筛选，能有助于获得准确的结果。筛选有2种方式，一种是通过添加参数来筛选。另一种是对 query 参数添加字段定义来筛选。

**参数过滤**：

常用的参数为日期参数

* -days 100：100天内的检索条件命中
* -maxdate：时间段最大日期，格式可以是 2001，2001/01/01
* -mindate：时间段最小日期，格式可以是 2001，2001/01/01
* -datetype: PDAT 发布时间；MDAT 修改时间

**字段定义**：

-query参数: esearch 和 efilter 通过添加该参数进行检索或过滤相应的关键词，query参数可以添加 filter，不同数据库可选择的filter可以用下面命令来查看。

```bash
# 查看sra数据库的filter
$ einfo -db sra -fields

# 查看sra数据库field，用xtract格式化输出
$ einfo -db sra | xtract -pattern Field -element Name Description
```

-query 检索词可以用逻辑关键词构建逻辑关系：

* AND：和
* OR：或
* NOT：非
* ()：共同定义

-query 检索词还可以使用正则表达式：

* "\*"星号表示任意个匹配的字符

**-query 字符串定义**

建议用 '...' 定义 -query 值，里面的各个关键词用 AND/OR 等建立逻辑关系。关键词组添加 filter 时用 "..." 定义。如：

-query '"salmonella paratyphi A"[ORGN] AND latest[filter] AND "complete genome[filter]"'

### 4. 具体示例

下面我们就用几个简单的例子来尝试 edirect 的功能。结合 Linux 其他命令行工具，可以极大的提升我们在命令行界面下检索和抓取 NCBI 数据库的效率。

**查看2014~2016年发表的霍乱弧菌CTX相关文献的摘要**

用关键词 vibrio cholerae CTX 检索 pubmed 数据库，将搜索到的文献的摘要保存到文件中，方便查看。由于是文本格式的内容，方便后期进一步处理，比如提取单词进行词频计算，用翻译API将摘要进行自动翻译等等。

```bash
# 示例为了减少命中数据，设置了发布日期的 -mindate 和 maxdata
$ esearch -db pubmed -mindate 2014 -maxdate 2016 -datetype PDAT \
> -query "vibrio cholerae CTX" | efetch -format abstract > abstract.txt

# 日期设置的另一种写法
$ esearch -db pubmed -query "vibrio cholerae CTX AND \
> (2014[pdat]:2016[pdat])" | efetch -format abstract > abstract.txt
```

**查看文章作者**

常常会对某一类文章，或者某一个期刊的文章作者进行关注，这里根据检索关键词，格式化输出文章的ID和作者列表。这也是 xtract 最常见的应用之一。

```bash
# 检索关键词，以最近300天为范围
$ esearch -db pubmed -query "listeria monocytogenes AND prfA" -days 300 | \
> efetch -format xml | xtract -pattern PubmedArticle \
> -element MedlineCitation/PMID -block Author -sep " " \
> -element LastName,Initials
```

**查看近300天以来发表的文章中，作者中有叫Jim的所有文章题目**

利用 pubmed 数据库的 filter [FULL], 查询所有作者列表里有叫Jim的文章标题，这里采用 medline 格式输出，如果要进一步做格式化可以用xml格式再用xtract重新排版打印。

```bash
$ esearch -db pubmed -days 300 -query "Jim[FULL]" | efetch -format medline | \
> grep TI | awk -F'- ' '{print $2}'
```

**查看2016年CDC发布的所有用miseq测序的沙门菌文库制备方法**

想要了解 CDC 开展的 miseq 沙门菌测序时采用的文库构建方法，可以利用 sra 数据库的 runinfo信息输出获得。

```bash
$ esearch -db sra -query "salmonella[OGRN] miseq[PLAT] EDBL-CDC" \
> -mindate 2016/01/01 -maxdate 2016/12/31 -datetype PDAT | \
> efetch -format runinfo | cut -d ',' -f 12 > library.txt
```

**绘制单增李斯特菌hlyA蛋白质进化树图**

```bash
$ esearch -db protein -query "listeriolysin hlyA NOT partial" |
> efetch -format fasta > hlyA.fasta
$ mafft hlyA.fasta > hlyA.mafft.fasta
$ raxml -m GTRGAMMA -p 12345 -s hlyA.mafft.fasta -f a -x 12345 -N 100 -T 40 -n hlyA
```

**对10年内发表的检索文献的作者排序**

```bash
$ esearch -db pubmed -days 3650 -datetype PDAT -query "hlyA" | \
> efetch -format xml | \
> xtract -pattern Author -element LastName | \
> sort-uniq-count-rank
```

**查看 taxonomy**

taxonomy 数据库的 filter [NXLV] 可以直接输出species级别

```bash
# 获得 taxonomy 数据库里 lis 或 sal 开头的分类名
$ esearch -db taxonomy -query "lis* OR sal*" | efetch -db taxonomy -format txt

# 获得 species 级的所有沙门菌
$ esearch -db taxonomy -query "Salmonella[NXLV]" | efetch -db taxonomy -format txt
```

**下载所有甲型副伤寒沙门菌基因组序列**

```bash
$ esearch -db nuccore -query "Salmonella[porgn:__txid54388] complete genome" | \
> efetch -db nuccore -format fasta > SPA.fasta
```

**获取SRA数据库中目的数据的测序样本地理信息**

```bash
$ esearch -db sra -query "China" -days 1 | elink -target biosample | \
> efetch -format docsum | xtract -pattern DocumentSummary -element Accession \
> -block Attribute -if Attribute@attribute_nam \
> -equals lat_lon -element Attribute
```

**获得某个bioproject下的所有sra run的信息**

```bash
$ eserach -db bioproject -query "PRJNA209644" | elink -target sra | \
> efetch -format docsum | xtract -pattern DocumentSummary -ACC @acc \
> -block DocumentSummary -element "&ACC"
```

**下载assembly数据库的refseq序列**

```bash
$ esearch -db assembly -query "Salmonella enterica[ORGN] AND 2010[GRLS]" | \
> efetch -format docsum | xtract -pattern DocumentSummary \
> -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}' | \
> xargs wget
```

**查看assembly meta信息**

NCBI 的 Assembly 数据库包含了 RefSeq 数据库的基因组数据，选择相应的物种可以获得该物种的基因组拼接序列，有完成图也有非完成图。

```bash
$ esearch -db assembly -query '"salmonella paratyphi A"[ORGN] AND latest[FILT]' | \
> efetch -format docsum | xtract -pattern Meta -block Stats \
> -tab '\n' -sep '\t' -element Stat@category Stat | sort -ur > meta.txt
```

```python
import pandas as pd
import seaborn as sns

df = pd.read_table("meta.txt")
sns.jointplot(x="contig_count", y="scaffold_n50", data=df)
```

![plot](../../assets/images/C02/03/plot.png)

**获取SRA菌株信息**

对某个物种开展进化研究时，我们往往会竟可能收集公共数据库中的菌株背景信息。下面这个例子搜索某个 bioproject，然后获得该研究的所有测序数据的 sra 背景信息。通过 \<SAMPLE_ATTRIBUTE\> 元素下的各个 TAG 和 VALUE 来提取所需信息。

```bash
$ esearch -db bioproject -query 273159 | elink -target sra | \
> efetch -format xml | xtract -pattern EXPERIMENT_PACKAGE -tab "," \
> -element RUN@accession -block SAMPLE_ATTRIBUTE -tab "\t" -sep ":" \
> -element TAG,VALUE > sample.txt
```

```bash
$ head sample.txt
```

这样输出的问题在于各个测序数据的 SAMPLE_ATTRIBUTE 的 TAG 顺序是乱序的，每个 SRR 输出的顺序不同导致各个字段要重新排序，比较麻烦。如果我们只需要其中一个字段，可以用 -if *** -eqauls *** 来过滤。比如下面代码我们输出菌株名称：

```bash
$ esearch -db bioproject -query 273159 | elink -target sra | \
> efetch -format xml | xtract -pattern EXPERIMENT_PACKAGE -tab "," \
> -element RUN@accession -block SAMPLE_ATTRIBUTE -if TAG -equals "strain" \
> -element VALUE | sort > strain-name.csv
```

`-if`在一行输出中不能反复使用，如果要同时对多个样品信息字段分别输出，我们采用变通的方法，先把 xml 数据保存到本地，然后循环分别输出各个所需内容，用 sort 排序，再用 paste 合并。

```bash
# 分别输出 csv 格式的数据
$ for i in serovar strain isolate_name_alias collection_date isolation_source; \
> do xtract -input bioproject.xml -pattern EXPERIMENT_PACKAGE -tab "," \
> -element RUN@accession -block SAMPLE_ATTRIBUTE -if TAG -equals $i -element \
> VALUE | sort > $i; done
# 合并 csv
$ paste -d , strain serovar isolate_name_alias collection_date > sample.csv
```

### Reference

[^1]: [Entrez Direct: E-utilities on the UNIX Command Line](http://www.ncbi.nlm.nih.gov/books/NBK179288)
[^4]: [Unix Pipeline Wikipedia](https://en.wikipedia.org/wiki/Pipeline_(Unix))
[^2]: [Medline Element Description](https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html)
[^3]: [XML](http://www.w3school.com.cn/xml/)

[NCBI]: https://www.ncbi.nlm.nih.gov/
[Entrez]: https://www.ncbi.nlm.nih.gov/gquery/
