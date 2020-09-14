# Mummer

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容介绍"
    Mummer 是一个历史悠久的快速细菌基因组比对工具，通过小片段`MUMs`比对到参考基因组，获得序列的SNP位点差异。因为其快速的特点，有许多细菌分析流程中都整合了Mummer。

Mummer算法是[Deltcher]()提出的一种基于后缀树数据结构的方法。大致原理是首先识别比对序列和参考序列的最大匹配片段（Maximal Unique Matches, MUMs），而该片段在比对序列和参考序列中是唯一的，不包含在其他更长MUMs中。对MUMs进行排序，然后对空位计算。Mummer适合对同种或近源基因组的序列比对，共线性分析，SNPs分析等。

## Mummer3.23

Mummer3.23是2003年发布的版本，一直沿用到今天。

### 安装

**conda安装**

```bash
# 新建虚拟环境
$ conda create -n mummer
# 进入虚拟环境
$ conda activate mummer
# 安装v3.23版本的mummer
(mummer)$ conda install mummer
# 安装gnuplot，v3版本的mummer只支持老版本的gnuplot
(mummer)$ conda install gnuplot=4.6.0
```

**源代码安装**

```bash
# 下载源代码
$ wget http://nchc.dl.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
# 解压缩
$ tar zxvf MUMmer3.23.tar.gz
$ cd MUMmer3.23
# 编译安装
$ make install
$ sudo ln -s `pwd`/bin/* /usr/local/sbin/
```

注意，编译需要的工具及版本。如果编译失败可能是缺乏依赖包或者版本太新：

- make (GNU make 3.79.1)
- perl (PERL     5.6.0)
- sh   (GNU sh   1.14.7)
- csh  (tcsh     6.10.00)
- g++  (GNU gcc  2.95.3)
- sed  (GNU sed  3.02)
- awk  (GNU awk  3.0.4)
- ar   (GNU ar   2.9.5)
- fig2dev (fig2dev 3.2.3)
- gnuplot (gnuplot 4.0)
- xfig    (xfig    3.2)

### 使用 Mummer3

本节将介绍执行mummer的基础知识，以及各个程序的用法及适用领域。应用程序详细信息，请参阅软件文档。

MUMmer软件包中最常用的五个程序是：

- mummer: 核心计算程序
- nucmer: 比对核酸的流程perl脚本
- promer: 比对氨基酸的流程perl脚本
- run-mummer1: shell分析流程脚本，
- run-mummer3: shell分析流程脚本，

其他功能程序

- show-coords: 将 .delta 数据文件生成 coordinates 数据输出
- show-snp
- show-aligns
- show-tiling
- show-diff
- mgaps: 将MUMs聚类
- annotate: 一个比对结果图片可视化的统一输出修正工具
- combineMUMs:
- delta-filter
- dnadiff
- exact-tandems
- mapview
- mummerplot
- repeat-match

**mummer**

mummer 是核心算法程序。其他的程序流程运行脚本。

参数：

- -mum: 计算唯一的最大匹配，不管是在参考序列还是比对序列中。
- -mumcand: 同 mumreference
- -mumreference: 计算比对到参考序列中的最大匹配，不考虑比对序列唯一性，这是默认设置
- -maxmatch: 计算所有最大匹配，不论他们是否是唯一性的
- -n: 设置匹配的字符，默认是 a/c/g/t，可以是大写和小写
- -l: 设置最小匹配序列长度，默认是20
- -b: 计算正义与翻译互补匹配
- -r: 只计算反义互补匹配
- -s: 显示匹配的信息
- -c: 报告反义互补匹配的比对序列位点信息
- -F: 强制使用4列输出格式
- -L: 在头部信息中显示比对序列的长度
- -h: 查看帮助信息

```bash
# 将query.fna比对到参考序列ref.fna
(mummer)$ mummer ref.fna query.fna
```


```bash
# 比较2个完成基因组
(mummer)$ nucmer
```

## Mummer4

!!! tip "注意事项"
    随着时代发展，新版的 Mummer4 尝试整合了一些新功能和多核性能的开放，但是目前还处于测试版，大部分的第三方分析流程也没有更新到 Mummer4。

### 安装

**conda安装**

```bash
$ conda create -n mummer4
$ conda activate mummer4
(mummer4)$ conda install mummer4
(mummer4)$ conda install gnuplot=4.6.0
```

**源代码安装**

```bash
$ wget https://github.com/gmarcais/mummer/releases/download/v4.0.0beta/mummer-4.0.0beta2.tar.gz
$ tar zxf mummer-4.0.0beta2.tar.gz && cd mummer-4.0.0beta2
$ ./configure
$ make && sudo make install
$ sudo cp .libs/libumdmummer.so.0.0.0 /lib/x86_64-linux-gnu/
$ sudo chmod 644 /lib/x86_64-linux-gnu/libumdmummer.so.0.0.0
$ sudo ln -s /lib/x86_64-linux-gnu/libumdmummer.so.0.0.0 /lib/x86_64-linux-gnu/libumdmummer.so.0
```

### 使用Mummer4

```bash
# 新增参数 -t 表示多线程计算
(mummer4)$ nucmer --maxmatch -t 40 -p test ref.fa query.fa
(mummer4)$ mummerplot -f --postscript --large test.delta
```

### Dot 工具做交互可视化

用 mummerplot 绘图是使用gnuplot绘制成共线性可视化图。mummer3 比较老，因此一些版本支持上会比较落后。这里推荐 [Dot](https://github.com/dnanexus/dot) 工具做交互式的可视化。[Dot](https://github.com/dnanexus/dot) 工具是 dnanexus 开发的一个在线可视化 nucmer 进行基因组共线性比较的工具，可以在线分析，也可以用浏览器在本地运行。

```bash
(mummer)$ git clone https://github.com/dnanexus/dot
(mummer)$ cd dot
(mummer)$ python DotPrep.py --delta compare.delta --out compare
(mummer)$ chromium index.html
```

上传 DotPrep.py 生成的 compare.coords 和 compare.coords.idx 文件到 https://dnanexus.github.io/dot 或者本地浏览器打开 index.html 页面，点击提交，获得可视化结果。

## 参考

1. [如何使用MUMmer比对大片段序列](https://vip.biotrainee.com/d/243-%E5%BA%94%E8%AF%A5%E6%98%AF%E6%9C%80%E8%AF%A6%E7%BB%86%E7%9A%84mummer%E4%B8%AD%E6%96%87%E4%BD%BF%E7%94%A8%E8%AF%B4%E6%98%8E)
