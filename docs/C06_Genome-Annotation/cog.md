# 基因功能注释

>### 同源基因概念

>**Homolog**: A gene related to a second gene by descent from a common ancestral DNA sequence. The term "homolog" may apply to the relationship between genes separated by speciation (ortholog), OR to the relationship betwen genes originating via genetic duplication (see paralog).

>**同源（基因）**：指遗传上来自于某一共同祖先DNA序列的基因。“同源”包含两种情况，一是物种间的同源（ortholog）；另一种是同一物种内由于基因复制、分离导致的同源（paralog）。

>Ortholog: Orthologs are genes in different species that have evolved from a common ancestral gene via speciation. Orthologs often (but certainly not always) retain the same function(s) in the course of evolution. Thus, functions may be lost or gained when comparing a pair of orthologs.

>**直系／直向／垂直同源（基因）**：指位于不同的物种间的，在物种形成过程中源自某一共同祖先的基因。从进化的角度来看，这类基因通常具有相同的功能，但并非绝对。所以，当我们比较一对直系同源基因时，可能会发现有的已经丧失了固有的功能或者进化出了新的功能。

>Paralog: Paralogs are genes produced via gene duplication within a genome. Paralogs typically evolve new functions or else eventually become pseudogenes.

>**旁系／横向／并行同源（基因）**：指在一个特定的基因组中由于基因复制产生的同源。通常比较这些基因时，它们可能已经彼此具有了新的（不同）功能，也可能已经成为假基因了。

## 1. COG

**COG**是**C**luster of **O**rthologous **G**roups of proteins（蛋白相邻类的聚簇）的缩写。构成每个COG的蛋白都是被假定为来自于一个祖先蛋白，并且因此或者是orthologs或者是paralogs。

- **Orthologs**是指来自于不同物种的由垂直家系（物种形成）进化而来的蛋白，并且典型的保留与原始蛋白有相同的功能。
- **Paralogs**是那些在一定物种中的来源于基因复制的蛋白，可能会进化出新的与原来有关的功能。

**COG**是通过把所有完整测序的基因组的编码蛋白一个一个的互相比较确定的。在考虑来自一个给定基因组的蛋白时，这种比较将给出每个其他基因组的一个最相似的蛋白。这些基因的每一个都轮番的被考虑。如果在这些蛋白（或子集）之间一个相互的最佳匹配关系被发现，那么那些相互的最佳匹配将形成一个COG。这样，一个COG中的成员将与这个COG中的其他成员比起被比较的基因组中的其他蛋白更相像，尽管如果绝对相似性比较的。最佳匹配原则的使用，没有了人为选择的统计切除的限制，这就兼顾了进化慢和进化快的蛋白。然而，还有一个加的限制就是一个COG必须包含来自于3个种系发生上远的基因组的一个蛋白。


### COG的分类


INFORMATION STORAGE AND PROCESSING

[J] Translation, ribosomal structure and biogenesis

[A] RNA processing and modification

[K] Transcription

[L] Replication, recombination and repair

[B] Chromatin structure and dynamics

CELLULAR PROCESSES AND SIGNALING

[D] Cell cycle control, cell division, chromosome partitioning

[Y] Nuclear structure

[V] Defense mechanisms

[T] Signal transduction mechanisms

[M] Cell wall/membrane/envelope biogenesis

[N] Cell motility

[Z] Cytoskeleton

[W] Extracellular structures

[U] Intracellular trafficking, secretion, and vesicular transport

[O] Posttranslational modification, protein turnover, chaperones

METABOLISM

[C] Energy production and conversion

[G] Carbohydrate transport and metabolism

[E] Amino acid transport and metabolism

[F] Nucleotide transport and metabolism

[H] Coenzyme transport and metabolism

[I] Lipid transport and metabolism

[P] Inorganic ion transport and metabolism

[Q] Secondary metabolites biosynthesis, transport and catabolism

POORLY CHARACTERIZED

[R] General function prediction only

[S] Function unknown

```bash
# 2003版 COGs
$ wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/myva
$ makeblastdb -in myva -dbtype prot -title cog -parse_seqids \
> -out PATH/TO/DB/cog -logfile PATH/TO/DB/cog.log

# 2014版 COGs
$
```

### 2014版 COGs的数据文件

>**genomes2003-2014.tab**
包含基因组信息的列表文件，格式采用：
`<genome-code> <ncbi-taxid> <ncbi-ftp-name>`
Acihos   	933801	  Acidianus_hospitalis_W1_uid66875

>**cognames2003-2014.tab**
COG注释信息的列表
`<COG-id> <functional-class> <COG-annotation>`
COG0001	     H	    Glutamate-1-semialdehyde aminotransferase

>**fun2003-2014.tab**
COG功能分类和描述
`<class-id> <description>`
J	Translation, ribosomal structure and biogenesis

>**cog2003-2014.csv**
同源基因域，以逗号分割
`<domain-id>, <genome-name>, <protein-id>,<protein-length>,<domain-start>, <domain-end>, <COG-id>, <membership-class>,`
333894695,Alteromonas_SN2_uid67349,333894695,427,1,427,COG0001,0,

>**prot2003-2014.tab**
所有蛋白质accession编号对应COGs的域名
`<protein-id> <refseq-acc>`
103485499	YP_615060

>**prot2003-2014.gi2gbk.tab**
所有蛋白质accession编号对应COGs的域名和Genbank编号
`<protein-id> <refseq-acc> <genbank-acc>`
103485499	YP_615060	ABF51727

>**prot2003-2014.fa.gz**
fasta格式的所有序列以及信息
`>gi|118430838|ref|NP_146899.2| putative mercury ion binding protein
[Aeropyrum pernix K1]
MIIFKRHSQAILFSHNKQEKALLGIEGMHCEGCAIAIETALKNVKGIIDTKVNYSRGSAI
VTFDDTLVSINDILEHYIFKVPSNYRAKLVSFIS`

## 2. KEGG

[KEGG](https://www.genome.jp/kegg/)是日本构建的为基因组、酶促途径以及生物化学物质在线数据库。它主要构建了代谢关系中的分子相互作用网络以及具体生物所有的变化形式。

系统信息
- PATHWAY — 细胞与生物机能的路径图
- MODULE — 基因模组与功能单位
- BRITE — 生物实体的阶层分类

基因组信息
- GENOME — 基因组全体
- GENES — 基因组全域的基因与蛋白列表
- ORTHOLOGY — 基因组全域的基因直系同源组别

化学信息
- COMPOUND, GLYCAN — 化合物与聚糖
- REACTION, RPAIR, RCLASS — 化学反应
- ENZYME — 酶命名法

健康信息
- DISEASE — 人类疾病
- DRUG — 已批准药物
- 环境 — 生药及与健康相关的物质

## 3. EggNOG

[EggNOG](http://eggnogdb.embl.de/#/app/home)

[EggNOG-mapper](https://github.com/jhcepas/eggnog-mapper)

### 3.1 本地安装 EggNOG-mapper

```bash
# 安装
$ conda create -n eggnog
$ conda activate eggnog
(eggnog)$ conda install eggnog-mapper
# 我们只关心细菌和病毒的基因组注释，因此只下载这2个数据库
# 各个单独数据库参见 http://beta-eggnogdb.embl.de/download/eggnog_4.5/hmmdb_levels
(eggnog)$ cd PATH/TO/DATA
(eggnog)$ download_eggnog_data.py -yq --data_dir `pwd` bact
```

### 3.2 注释

emapper 的默认算法是hmm，如果想用diamond来做cluster，可以用`-m diamond`参数

```bash
# HMM 注释
# 针对细菌基因组数据库，使用40个CPU内核，数据库运行在内存中
(eggnog)$ emapper.py -i mygenome.fa -d bact --data_dir PATH/TO/DATA \
> --cpus 40 --usemem -o mygenome_bact_result
# 针对病毒基因组数据库
(eggnog)$ emapper.py -i mygenome.fa -d viruses --data_dir PATH/TO/DATA -o mygenome_virus_result
# Diamond 注释
(eggnog)$ emapper.py -i mygenome.fa -d bact -m diamond --data_dir PATH/TO/DATA -o mygenome_bact_result
```

### 3.3 结果

- \*.emapper.hmm_hits
- \*.emapper.seed_orthologs
- \*.emapper.annotations

hmm_hits 文件

## 参考资料

1. [NCBI的COG介绍](https://www.plob.org/article/1000.html)
2. [Expanded microbial genome coverage and improved protein family annotation in the COG database](https://www.ncbi.nlm.nih.gov/pubmed/25428365)
3. Fitch WM (1970) "Distinguishing homologous from analogous proteins". Syst Zool, 19:99-113.
4. Li WH, Yang J, Gu X (2005) "Expression divergence between duplicate genes". Trends Genet, 21:602-607.
5. http://blog.sciencenet.cn/blog-217859-280960.html
6. http://homepage.usask.ca/~ctl271/857/def_homolog.shtml
