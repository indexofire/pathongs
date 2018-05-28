# Metagenomics数据分析

病原微生物`Metagenomics`根据研究内容，可以分为几个方向：`whole genomes`测序，例如`shotgun metagenomics`寻找所有序列;`16s rRNA`测序，例如研究群菌结构；`Transcriptome`测序，研究整个样品中转录组。

**Metagenomics 数据分析的基本流程：**

1. Trimming: 切除引物序列，barcode序列等
2. Binnings:
3. OTUs:
4. Chimeras

**常用工具：**

1. MEGAN
2. mothur
3. QIIME
4. AXIOME & CloVR

---

## 1. Meta-BEETL

`BEETL` 是一套做 `Burrows-Wheeler Transform` 的工具集，[Meta-BEETL] 是 `BEETL` 下的一个子项目，用来分析基于 Shotgun metagenomics 的微生物数据。

### 1. 已整理的数据库

作者已经将研究所用到的数据库整理后提交到 Amazon S3 服务器上，可以直接下载。总数据量在24G左右。

```bash
~$ mkdir BeetlMetagenomeDatabase && cd BeetlMetagenomeDatabase
~/BeetlMetagenomDatabase$ for i in `curl https://s3.amazonaws.com/metaBEETL | grep -oP "[^>]*bz2"` ; \
> do wget https://s3.amazonaws.com/metaBEETL/$i & done
```

静态文件的 URL 地址，也可以直接下载。

```
https://s3.amazonaws.com/metaBEETL/filecounter.csv.bz2
https://s3.amazonaws.com/metaBEETL/headerFile.csv.bz2
https://s3.amazonaws.com/metaBEETL/names.dmp.bz2
https://s3.amazonaws.com/metaBEETL/ncbiFileNumToTaxTree.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B00.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B01.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B02.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B03.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B04.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B05.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-B06.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C00.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C01.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C02.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C03.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C04.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C05.bz2
https://s3.amazonaws.com/metaBEETL/ncbiMicros-C06.bz2
```
下载完成后解压缩。

```
~/BeetlMetagenomDatabase$ bunzip2 *.bz2
~/BeetlMetagenomDatabase$ export METAGENOME_DATABASE_PATH=`pwd`
```

### 2. 自己建立数据库

处于实际目的，有时候我们不需要这么大的数据库，或者需要更多其他的数据加入到数据库，那就需要自行建立数据库。

#### 2.1 安装 SeqAn

[SeqAn](www.seqan.de) 是一个高效的 C++ 算法和数据结构库，用于分析生物序列等用途。访问 [SeqAn](www.seqan.de) 下载并安装。安装完成后编译BEETL时可以添加`--with-seqan`参数。

```bash
# 安装依赖包
~$ sudo apt-get install g++ cmake python zlib1g-dev libbz2-dev libboost-dev
~$ cd ~/apps
~/apps$ git clone https://github.com/seqan/seqan.git && cd seqan
~/apps/seqan$ git checkout -b develop origin/develop
```

#### 2.2 下载目标序列

下载需要的序列，比如NCBI上细菌的基因组 ref 数据。

```bash
~/BeetlMetagenomDatabase$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz -P downloads
~/BeetlMetagenomDatabase$ mkdir all.fna
~/BeetlMetagenomDatabase$ tar xzf downloads/all.fna.tar.gz -C all.fna/
```

#### 2.3 创建去除质粒序列的单基因组序列

```bash
~/BeetlMetagenomDatabase$ metabeetl-db-genomesToSingleSeq -f all.fna/*/*.fna -s singleSeqGenomes
~/BeetlMetagenomDatabase$ cd singleSeqGenomes
```

#### 2.4 创建BWT

修改`metabeetl-db-arrayBWT.sh`文件里数据库路径，创建G\_\* 和 G\_\*\_rev 文件的BWTs，用qsub来向系统提交任务。也可以运行`metabeetl-db-makeBWTSkew`来一个个创建。

```bash
~/BeetlMetagenomDatabase$ cp `which metabeetl-db-arrayBWT.sh` .
~/BeetlMetagenomDatabase$ vim metabeetl-db-arrayBWT.sh
# 如果文件是 G_1 - G_500，那么n=1-500
~/BeetlMetagenomDatabase$ qsub -t n metabeetl-db-arrayBWT.sh

# 或者
~/BeetlMetagenomDatabase$ (echo -n "all: " ; for i in G_*; do echo -n " bwt_${i}-B00"; done ; echo -e "\n" ; \
> for i in G_*; do echo "bwt_${i}-B00: ${i}"; echo -e "\tmetabeetl-db-makeBWTSkew ${i} ${i}\n" ; done ) > Makefile
~$ make -j
```

#### 2.5 合并序列

在一台内存>60G的机器上将序列加载到内存中，并吧所有文件合并。

```bash
~/BeetlMetagenomDatabase$ for pileNum in `seq 0 5`; do metabeetl-db-mergeBacteria $pileNum ncbiMicros <( ls G_* ) ; done
```
生成的文件名称类似：ncbiMicros-A0\*, -B0\* and -C0\*  

- -A0\* 包含每个BWT的位置
- -B0\* BWTs
- -C0\* 包含每个BWT位置的文件号

   文件 -B0\* and -C0\* 用于计数算法

#### 2.6 转换BWT文件格式

可选操作：将BWT文件转换成RLE BWT格式，运行可以更快。

```bash
~/BeetlMetagenomDatabase$ for pileNum in `seq 0 5`; do \
> mv ncbiMicros-B0${pileNum} ncbiMicros-B0${pileNum}.ascii ;  \
> beetl-convert \  
> --input-format=bwt_ascii \  
> --output-format=bwt_rle \  
> -i ncbiMicros-B0${pileNum}.ascii \  
> -o ncbiMicros-B0${pileNum} ; \  
> done
```

#### 2.7 下载NCBI Taxoonmy

```bash
~$ cd tmp
~/tmp$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
~/tmp$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
~/tmp$ tar xzf taxdump.tar.gz names.dmp nodes.dmp
~/tmp$ gunzip gi_taxid_nucl.dmp.gz
```

用`metabeetl-db-findTaxa`脚本查找数据库中的taxonomic树。
   Use the metabeetl-db-findTaxa script to find the taxonomic tree corresponding to the file numbers in the database.  
   You will need the headerFile produced by running "metabeetl-db-genomesToSingleSeq" and
   fileCounter created during the merging of the bacterial reference genomes.  
   Finally, you get for each file number in the database a taxonomic tree with the taxonomic ids.  
   There will be some 0 in the taxonomic tree. This is a taxonomic id which could not be
   matched to: Superkingdom, Phylum, Order, Class, Family, Genus, Species or Strain.  
   Sometimes there are just missing taxa in the taxonomy. We supplement this with the file `metaBeetlExtraNames.dmp` below.

```bash
~/BeetlMetagenomDatabase$ metabeetl-db-findTaxa \
> -nA downloads/names.dmp \
> -nO downloads/nodes.dmp \
> -nG downloads/gi_taxid_nucl.dmp \
> -h singleSeqGenomes/headerFile.csv \
> -f singleSeqGenomes/filecounter.csv > ncbiFileNumToTaxTree

~/BeetlMEtagenomDatabase$ ( grep scientific downloads/names.dmp ; cat ${BEETL_INSTALL_DIR}/share/beetl/metaBeetlExtraNames.dmp ) > metaBeetlTaxonomyNames.dmp
```

#### 2.8 计算 normalisation 因子

```bash
~/BeetlMetagenomDatabase$ mkdir normalisation
~/BeetlMetagenomDatabase$ cd normalisation
~/BeetlMetagenomDatabase/normalisation$ touch normalize.sh
```

```bash
# normalize.sh
for genome in ../singleSeqGenomes/G_*; do
(
    genomeNum=`basename ${genome}`
    mkdir ${genomeNum}
    cd ${genomeNum}
    for i in ../../singleSeqGenomes/bwt_${genomeNum}-B0?; do ln -s ${i} ; done
    for i in `seq 0 6`; do touch bwt_${genomeNum}-B0${i}; done
    time beetl-compare --mode=metagenomics -a bwt_${genomeNum} -b ../../singleSeqGenomes/ncbiMicros -t ../../ncbiFileNumToTaxTree -w 20 -n 1 -k 50 --no-comparison-skip -v &> out.step1
    rm -f BeetlCompareOutput/cycle51.subset*
        touch empty.txt
    time cat BeetlCompareOutput/cycle*.subset* | metabeetl-convertMetagenomicRangesToTaxa ../../ncbiFileNumToTaxTree ../../singleSeqGenomes/ncbiMicros ../../metaBeetlTaxonomyNames.dmp empty.txt 20 50 - &> out.step2
    cd ..
) &
done ; wait

for i in `seq 1 2977`; do echo "G_$i"; X=`grep -P "G_${i}$" ../singleSeqGenomes/filecounter.csv |cut -f 1 -d ','`; TAX=`grep -P "^${X} " ../ncbiFileNumToTaxTree | tr -d "\n\r"` ; echo "TAX(${X}): ${TAX}."; TAXID=`echo "${TAX}" | sed 's/\( 0\)*$//g' |awk '{print $NF}'`; echo "TAXID=${TAXID}"; COUNTS=`grep Counts G_${i}/out.step2 | head -1`; echo "COUNTS=${COUNTS}"; MAIN_COUNT=`echo "${COUNTS}  " | sed "s/^.* ${TAXID}:\([0-9]*\) .*$/\1/ ; s/Counts.*/0/"` ; echo "MAIN_COUNT=${MAIN_COUNT}" ; SUM=`echo "${COUNTS}  " | tr ' ' '\n' | sed 's/.*://' | awk 'BEGIN { sum=0 } { sum+=$1 } END { print sum }'` ; echo "SUM=$SUM"; PERCENT=`echo -e "scale=5\n100*${MAIN_COUNT}/${SUM}" | bc` ; echo "PERCENT=${PERCENT}" ; echo "FINAL: G_${i} ${TAXID} ${MAIN_COUNT} ${SUM} ${COUNTS}" ; done > r2977

grep FINAL r2977 > ../normalisation.txt
```

#### 2.9 创建数据库文件夹

```bash
~$ mkdir metaBeetlNcbiDb
~$ cd ~/metaBeetlNcbiDb
~$ ln -s ../ncbiFileNumToTaxTree
~$ ln -s ../normalisation.txt
~$ ln -s ../downloads/metaBeetlTaxonomyNames.dmp
~$ ln -s ../singleSeqGenomes/filecounter.csv
~$ ln -s ../singleSeqGenomes/headerFile.csv
~$ for i in ../singleSeqGenomes/ncbiMicros-[BC]0[0-5]; do ln -s $i ; done
```

### 3. 分析数据

#### 3.1 获取数据

用 SRS013948 这个人类肠道细菌组研究项目的数据作为例子。首先下载数据:

```bash
# http download
~/data$ wget http://downloads.hmpdacc.org/data/Illumina/throat/SRS013948.tar.bz2

# ftp download
~/data$ wget ftp://public-ftp.hmpdacc.org/Illumina/throat/SRS013948.tar.bz2

# 解压缩
~/data$ tar xjf SRS013948.tar.bz2
~/data$ ls -l
```

解压缩后可以看到下面3个文件：

```
SRS013948.denovo_duplicates_marked.trimmed.1.fastq
SRS013948.denovo_duplicates_marked.trimmed.2.fastq
SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq
```

#### 3.2 转换数据和合并数据

```bash
~/data$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.1.fastq \
> -o paddedSeq1.seq \
> --sequence-length=100
~/data$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.2.fastq \
> -o paddedSeq2.seq \
> --sequence-length=100
~/data$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq \
> -o paddedSeqSingleton.seq \
> --sequence-length=100
~/data$ cat paddedSeq1.seq paddedSeq2.seq paddedSeqSingleton.seq > SRS013948.seq
```

#### 3.3 以 metagenomic 模式运行 BEETL

创建 BWT 变换，然后运行 metagenomic 模式的 BEETL

```bash
~/data$ beetl-bwt -i SRS013948.seq -o bwt_SRS013948
~/data$ beetl-compare \
> --mode=metagenomics \
> -a bwt_SRS013948 \
> -b ${METAGENOME_DATABASE_PATH}/ncbiMicros \
> -t ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
> -w 20 \
> -n 1 \
> -k 50 \
> --no-comparison-skip
```

`beetl-compare` 命令会在文件夹 `BeetlCompareOutput` 里生成许多文件。k值设置越大，可以获得越多的信息数据，但是输出文件也会变得越大。

#### 3.4 图形化显示结果

`metabeetl-convertMetagenomicRangesToTaxa` 工具根据kmer匹配，从而获得基因组以及上一级分类的ID号，并生成文字与图形结果（采用Krona js可视化库来形成网页格式的报告文件）。

由于算法原因需要频繁读取BWT位置，因此数据库文件 `ncbiMicros-C0*` 最好放在读取速度比较块的磁盘扇区或者磁盘里如SSD硬盘。另外如果内存足够，可以用 `metabeetl-convertMetagenomicRangesToTaxa_withMmap` 工具，将这些比对文件读入内存中加快速度。

```bash
~/data$ cat BeetlCompareOutput/cycle*.subset* | \
> metabeetl-convertMetagenomicRangesToTaxa \
> ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
> ${METAGENOME_DATABASE_PATH}/ncbiMicros \
> ${METAGENOME_DATABASE_PATH}/metaBeetlTaxonomyNames.dmp \
> ${METAGENOME_DATABASE_PATH}/normalisation.txt \
> 20 50 - > metaBeetl.log
```

**生成的TSV（制表符分隔）文件各列含义：**

| 列 | 字段含义 |
| -- | -- |
| column 1 | Taxonomy Id |
| column 2 | Taxonomy level |
| column 3 | k-mer count |
| column 4 | k-mer count including children |

**TSV文件用途：**

- metaBeetl.tsv: raw k-mer counts for every leaves and ancestors.
- metaBeetl_normalised.tsv: some counts from ancestors are moved towards leaf items. The proportion thereof is pre-computed by aligning each individual genome to the full database.
- metaBeetl_normalised2.tsv: Only leaves
 of the taxonomy tree are kept, and counts are normalised relatively to genome sizes.

**Krona JS 结果文件**
- metaBeetl_krona.html
- metaBeetl_krona_normalised.html
- metaBeetl_krona_normalised2.html

## Reference

* http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3622627/
* https://github.com/BEETL/BEETL.git

[metaBEETL]: https://github.com/BEETL/BEETL


## 2. SUPRI

SUPRI 是主要用来对临床样品的 shotgun metagenomics 测序数据中寻找并鉴定病原的工具。它有2种模式：Fast mode 和 Comprehensive mode。前者可以快速对测序的 reads 进行 Mapping，找到细菌或病毒 reads。后者分析更为详细，除了细菌病毒外还比对了真菌，寄生虫等其他物种数据库，并且进行de novo assembly，对contig也进行 Mapping，虽然耗时更长，但可以获得覆盖度，taxonomic 分类等信息。

### 安装依赖工具

* [fastQValidator](http://genome.sph.umich.edu/wiki/FastQValidator)
* [Minimo (v1.6)](http://sourceforge.net/projects/amos/files/amos/3.1.0/)
* [Abyss (v1.3.5)](http://www.bcgsc.ca/platform/bioinfo/software/abyss)
* [RAPSearch (v2.12)](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
* [seqtk (v 1.0r31)](https://github.com/lh3/seqtk)
* [SNAP (v0.15)](http://snap.cs.berkeley.edu)
* [gt (v1.5.1)](http://genometools.org/index.html)
* [fastq](https://github.com/brentp/bio-playground/tree/master/reads-utils)
* [fqextract](https://gist.github.com/drio/1168330)
* [cutadapt (v1.2.1)](https://code.google.com/p/cutadapt/)
* [prinseq](http://prinseq.sourceforge.net)
* [dropcache](http://stackoverflow.com/questions/13646925/allowing-a-non-root-user-to-drop-cache)

### 安装 SUPRI

```bash
$ wget https://github.com/chiulab/surpi/releases/download/v1.0.18/surpi-1.0.18.tar.gz
$ tar zxf surpi-1.0.18.tar.gz -C surpi
```

### 建立数据库

```bash
$ mkdir -p SNAP_db
$ ~/apps/surpi/create_taxonomy_db.sh
```

### 运行SURPI

```bash
$ ~/apps/surpi/SURPI.sh -z input.fastq
$ ~/apps/surpi/go_input &
```

### 输出结果

当以 Comprehensive mode 运行后，会产生以下文件夹：

```
DATASETS_input
deNovoASSEMBLY_input
LOG_input
OUTPUT_input
TRASH_input
```

结果文件在`OUTPUT_input`文件夹里，“.annotated” 文件是 SAM 格式或者 -m8 格式的比对结果, taxonomic 信息也在文件中的最后一行中。 “.counttable” 文件是 tab 分割的总结文件。每一行表示不同层面的 taxonomic 注释信息。

## Reference

1. [SUPRI official](http://chiulab.ucsf.edu/surpi)


## Reference

1. [Introduction to Metagenomics Data Analysis](http://www.slideshare.net/ueb52/2013-08mg-training)
