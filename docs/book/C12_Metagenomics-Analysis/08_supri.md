# Supri

---

!!! Abstract "内容简介"
    Supri 是一个老牌的宏基因组序列扫描软件，创建它之初是用来做临床上病原微生物序列检索的。

SUPRI 是主要用来对临床样品的 shotgun metagenomics 测序数据中寻找并鉴定病原的工具。它有2种模式：Fast mode 和 Comprehensive mode。前者可以快速对测序的 reads 进行 Mapping，找到细菌或病毒 reads。后者分析更为详细，除了细菌病毒外还比对了真菌，寄生虫等其他物种数据库，并且进行de novo assembly，对contig也进行 Mapping，虽然耗时更长，但可以获得覆盖度，taxonomic 分类等信息。

## 1.安装依赖工具

先要安装依赖工具：

- [fastQValidator](http://genome.sph.umich.edu/wiki/FastQValidator)
- [Minimo v1.6](http://sourceforge.net/projects/amos/files/amos/3.1.0/)
- [Abyss v1.3.5](http://www.bcgsc.ca/platform/bioinfo/software/abyss)
- [RAPSearch v2.12](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
- [seqtk v 1.0r31](https://github.com/lh3/seqtk)
- [SNAP v0.15](http://snap.cs.berkeley.edu)
- [gt v1.5.1](http://genometools.org/index.html)
- [fastq](https://github.com/brentp/bio-playground/tree/master/reads-utils)
- [fqextract](https://gist.github.com/drio/1168330)
- [cutadapt v1.2.1](https://code.google.com/p/cutadapt/)
- [prinseq-lite.pl](http://prinseq.sourceforge.net)
- [dropcache](http://stackoverflow.com/questions/13646925/allowing-a-non-root-user-to-drop-cache)

## 2.安装 SUPRI

```bash
$ wget https://github.com/chiulab/surpi/releases/download/v1.0.18/surpi-1.0.18.tar.gz
$ tar zxf surpi-1.0.18.tar.gz -C surpi
```

## 3.建立数据库

```bash
$ mkdir -p SNAP_db
$ create_taxonomy_db.sh
```

## 4.运行SURPI

```bash
$ SURPI.sh -z input.fastq
$ ./go_input &
```

## 5.输出结果

当以 Comprehensive mode 运行后，会产生以下文件夹：

- DATASETS_input
- deNovoASSEMBLY_input
- LOG_input
- OUTPUT_input
- TRASH_input

结果文件在OUTPUT_input文件夹里，“.annotated” 文件是 SAM 格式或者 -m8 格式的比对结果, taxonomic 信息也在文件中的最后一行中。 “.counttable” 文件是 tab 分割的总结文件。每一行表示不同层面的 taxonomic 注释信息。

## 参考资料

1. [SUPRI official](https://chiulab.ucsf.edu/surpi/)
