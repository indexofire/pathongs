# Meta-BEETL

---

!!! Abstract "内容简介"
    `BEETL`是一套做`Burrows-Wheeler Transform`的工具集，`Meta-BEETL <#meta-beetl>`是`BEETL`下的一个子项目，用来分析基于 Shotgun metagenomics 的微生物数据。

## Meta-BEETL Metagenomics

### 下载数据库

#### 1. 已整理的数据库

作者已经将研究所用到的数据库整理后提交到 Amazon S3 服务器上，可以直接下载。总数据量在24G左右。

```bash
$ mkdir BeetlMetagenomeDatabase && cd BeetlMetagenomeDatabase
$ for i in `curl https://s3.amazonaws.com/metaBEETL | grep -oP "[^>]*bz2"` ; \
    > do wget https://s3.amazonaws.com/metaBEETL/$i & done
```

静态文件的 URL 地址，也可以直接下载。

- https://s3.amazonaws.com/metaBEETL/filecounter.csv.bz2
- https://s3.amazonaws.com/metaBEETL/headerFile.csv.bz2
- https://s3.amazonaws.com/metaBEETL/names.dmp.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiFileNumToTaxTree.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B00.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B01.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B02.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B03.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B04.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B05.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-B06.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C00.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C01.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C02.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C03.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C04.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C05.bz2
- https://s3.amazonaws.com/metaBEETL/ncbiMicros-C06.bz2

下载完成后解压缩。

```bash
$ bunzip2 *.bz2
$ export METAGENOME_DATABASE_PATH=`pwd`
```

#### 2. 自己建立数据库

处于实际目的，有时候我们不需要这么大的数据库，或者需要更多其他的数据加入到数据库，那就需要自行建立数据库。

**2.1 安装 SeqAn**

`SeqAn <www.seqan.de>`__ 是一个高效的 C++
算法和数据结构库，用于分析生物序列等用途。访问 `SeqAn <www.seqan.de>`__
下载并安装。安装完成后编译BEETL时可以添加\ ``--with-seqan``\ 参数。

```bash
# 安装依赖包
$ sudo apt-get install g++ cmake python zlib1g-dev libbz2-dev libboost-dev
$ cd ~/apps
$ git clone https://github.com/seqan/seqan.git && cd seqan
$ git checkout -b develop origin/develop
```

**2.2 下载目标序列**

下载需要的序列，比如NCBI上细菌的基因组 ref 数据。

```bash
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz -P downloads
$ mkdir all.fna
$ tar xzf downloads/all.fna.tar.gz -C all.fna/
```

**2.3 创建去除质粒序列的单基因组序列**

```bash
$ metabeetl-db-genomesToSingleSeq -f all.fna/*/*.fna -s singleSeqGenomes
$ cd singleSeqGenomes
```

**2.4 创建BWT**

修改\ ``metabeetl-db-arrayBWT.sh``\ 文件里数据库路径，创建G\_\* 和
G\_\*\_rev 文件的BWTs，用qsub来向系统提交任务。也可以运行\ ``metabeetl-db-makeBWTSkew``\ 来一个个创建。

```bash
$ cp `which metabeetl-db-arrayBWT.sh` .
$ vim metabeetl-db-arrayBWT.sh
$ metabeetl-db-arrayBWT.sh

# 或者使用循环
$ (echo -n "all: " ; for i in G_*; do echo -n " bwt_${i}-B00"; done ; echo -e "\n" ; \
> for i in G_*; do echo "bwt_${i}-B00: ${i}"; echo -e "\tmetabeetl-db-makeBWTSkew ${i} ${i}\n" ; done ) > Makefile
$ make -j
```

**2.5 合并序列**

在一台内存>60G的机器上将序列加载到内存中，并吧所有文件合并。

```bash
$ for pileNum in `seq 0 5`; do metabeetl-db-mergeBacteria $pileNum ncbiMicros <( ls G_* ) ; done
```

生成的文件名称类似：ncbiMicros-A0\*, -B0\* and -C0\*

- A0\* 包含每个BWT的位置
- B0\* BWTs
- C0\* 包含每个BWT位置的文件号

文件 -B0\* and -C0\* 用于计数算法

**2.6 转换BWT文件格式**

可选操作：将BWT文件转换成RLE BWT格式，运行可以更快。

```bash
$ for pileNum in `seq 0 5`; do \
> mv ncbiMicros-B0${pileNum} ncbiMicros-B0${pileNum}.ascii ;  \
> beetl-convert \  
> --input-format=bwt_ascii \  
> --output-format=bwt_rle \  
> -i ncbiMicros-B0${pileNum}.ascii \  
> -o ncbiMicros-B0${pileNum} ; \  
> done
```

**2.7 下载NCBI Taxoonmy**

```bash
$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
$ tar xzf taxdump.tar.gz names.dmp nodes.dmp
$ gunzip gi_taxid_nucl.dmp.gz
```

!!! infomation
    用`metabeetl-db-findTaxa`脚本查找数据库中的taxonomic树。
    You will need the headerFile produced by running "metabeetl-db-genomesToSingleSeq" and fileCounter created during the merging of the bacterial reference genomes.

|  Finally, you get for each file number in the database a taxonomic tree with the taxonomic ids.

|  There will be some 0 in the taxonomic tree. This is a taxonomic id which could not be matched to: Superkingdom, Phylum, Order, Class, Family, Genus, Species or Strain.

|  Sometimes there are just missing taxa in the taxonomy. We supplement this with the file ``metaBeetlExtraNames.dmp`` below.

```bash
$ metabeetl-db-findTaxa \
> -nA downloads/names.dmp \
> -nO downloads/nodes.dmp \
> -nG downloads/gi_taxid_nucl.dmp \
> -h singleSeqGenomes/headerFile.csv \
> -f singleSeqGenomes/filecounter.csv > ncbiFileNumToTaxTree

$ ( grep scientific downloads/names.dmp ; cat ${BEETL_INSTALL_DIR}/share/beetl/metaBeetlExtraNames.dmp ) > metaBeetlTaxonomyNames.dmp
```

**2.8 计算 normalisation 因子**

```bash
$ mkdir normalisation
$ cd normalisation
$ touch normalize.sh
$ vim normalize.sh
```

编辑脚本内容如下

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

for i in `seq 1 2977`; do echo "G_$i"; X=`grep -P "G_${i}$" ../singleSeqGenomes/filecounter.csv |cut -f 1 -d ','`; TAX=`grep -P "^${X} " ../ncbiFileNumToTaxTree | tr -d "\n\r"` ; echo "TAX(${X}): ${TAX}."; TAXID=`echo "${TAX}" | sed 's/\( 0\)*$//g' | awk '{print $NF}'`; echo "TAXID=${TAXID}"; COUNTS=`grep Counts G_${i}/out.step2 | head -1`; echo "COUNTS=${COUNTS}"; MAIN_COUNT=`echo "${COUNTS}  " | sed "s/^.* ${TAXID}:\([0-9]*\) .*$/\1/ ; s/Counts.*/0/"` ; echo "MAIN_COUNT=${MAIN_COUNT}" ; SUM=`echo "${COUNTS}  " | tr ' ' '\n' | sed 's/.*://' | awk 'BEGIN { sum=0 } { sum+=$1 } END { print sum }'` ; echo "SUM=$SUM"; PERCENT=`echo -e "scale=5\n100*${MAIN_COUNT}/${SUM}" | bc` ; echo "PERCENT=${PERCENT}" ; echo "FINAL: G_${i} ${TAXID} ${MAIN_COUNT} ${SUM} ${COUNTS}" ; done > r2977

grep FINAL r2977 > ../normalisation.txt
```

**2.9 创建数据库文件夹**

```bash
$ mkdir metaBeetlNcbiDb
$ cd ~/metaBeetlNcbiDb
$ ln -s ../ncbiFileNumToTaxTree
$ ln -s ../normalisation.txt
$ ln -s ../downloads/metaBeetlTaxonomyNames.dmp
$ ln -s ../singleSeqGenomes/filecounter.csv
$ ln -s ../singleSeqGenomes/headerFile.csv
$ for i in ../singleSeqGenomes/ncbiMicros-[BC]0[0-5]; do ln -s $i ; done
```

#### 3. 分析数据

**3.1 获取数据**

用 SRS013948 这个人类肠道细菌组研究项目的数据作为例子。首先下载数据:

```bash
# http download
$ wget http://downloads.hmpdacc.org/data/Illumina/throat/SRS013948.tar.bz2

# ftp download
$ wget ftp://public-ftp.hmpdacc.org/Illumina/throat/SRS013948.tar.bz2

# 解压缩
$ tar xjf SRS013948.tar.bz2
$ ls -l
```

解压缩后可以看到下面3个文件：

- SRS013948.denovo_duplicates_marked.trimmed.1.fastq
- SRS013948.denovo_duplicates_marked.trimmed.2.fastq
- SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq

**3.2 转换数据和合并数据**

```bash
$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.1.fastq \
> -o paddedSeq1.seq \
> --sequence-length=100

$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.2.fastq \
> -o paddedSeq2.seq \
> --sequence-length=100

$ beetl-convert \
> -i SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq \
> -o paddedSeqSingleton.seq \
> --sequence-length=100

$ cat paddedSeq1.seq paddedSeq2.seq paddedSeqSingleton.seq > SRS013948.seq
```

**3.3 以 metagenomic 模式运行 BEETL**

创建 BWT 变换，然后运行 metagenomic 模式的 BEETL

```bash
$ beetl-bwt -i SRS013948.seq -o bwt_SRS013948
$ beetl-compare \
> --mode=metagenomics \
> -a bwt_SRS013948 \
> -b ${METAGENOME_DATABASE_PATH}/ncbiMicros \
> -t ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
> -w 20 \
> -n 1 \
> -k 50 \
> --no-comparison-skip
```

`beetl-compare` 命令会在文件夹 `BeetlCompareOutput`里生成许多文件。k值设置越大，可以获得越多的信息数据，但是输出文件也会变得越大。

**3.4 图形化显示结果**


`metabeetl-convertMetagenomicRangesToTaxa`: 工具根据kmer匹配，从而获得基因组以及上一级分类的ID号，并生成文字与图形结果（采用Krona js可视化库来形成网页格式的报告文件）。由于算法原因需要频繁读取BWT位置，因此数据库文件 `ncbiMicros-C0*`最好放在读取速度比较块的磁盘扇区或者磁盘里如SSD硬盘。另外如果内存足够，可以用`metabeetl-convertMetagenomicRangesToTaxa_withMmap`工具，将这些比对文件读入内存中加快速度。

```bash
$ cat BeetlCompareOutput/cycle*.subset* | \
> metabeetl-convertMetagenomicRangesToTaxa \
> ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
> ${METAGENOME_DATABASE_PATH}/ncbiMicros \
> ${METAGENOME_DATABASE_PATH}/metaBeetlTaxonomyNames.dmp \
> ${METAGENOME_DATABASE_PATH}/normalisation.txt \
> 20 50 - > metaBeetl.log
```

**生成的TSV（制表符分隔）文件各列含义：**


| 列         | 字段含义                         |
|------------|----------------------------------|
| column 1   | Taxonomy Id                      |
| column 2   | Taxonomy level                   |
| column 3   | k-mer count                      |
| column 4   | k-mer count including children   |

**TSV文件用途：**

-  metaBeetl.tsv: raw k-mer counts for every leaves and ancestors.
-  metaBeetl\_normalised.tsv: some counts from ancestors are moved
   towards leaf items. The proportion thereof is pre-computed by
   aligning each individual genome to the full database.
-  metaBeetl\_normalised2.tsv: Only leaves of the taxonomy tree are
   kept, and counts are normalised relatively to genome sizes.

**Krona JS 结果文件** - metaBeetl\_krona.html -
metaBeetl\_krona\_normalised.html - metaBeetl\_krona\_normalised2.html

## 参考资料

- http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3622627/
- https://github.com/BEETL/BEETL.git
