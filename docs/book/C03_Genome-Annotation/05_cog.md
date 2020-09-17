# COG 数据库

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

>### 同源基因概念

>**Homolog**: A gene related to a second gene by descent from a common ancestral DNA sequence. The term "homolog" may apply to the relationship between genes separated by speciation (ortholog), OR to the relationship betwen genes originating via genetic duplication (see paralog).

>**同源（基因）**：指遗传上来自于某一共同祖先DNA序列的基因。“同源”包含两种情况，一是物种间的同源（ortholog）；另一种是同一物种内由于基因复制、分离导致的同源（paralog）。

>Ortholog: Orthologs are genes in different species that have evolved from a common ancestral gene via speciation. Orthologs often (but certainly not always) retain the same function(s) in the course of evolution. Thus, functions may be lost or gained when comparing a pair of orthologs.

>**直系／直向／垂直同源（基因）**：指位于不同的物种间的，在物种形成过程中源自某一共同祖先的基因。从进化的角度来看，这类基因通常具有相同的功能，但并非绝对。所以，当我们比较一对直系同源基因时，可能会发现有的已经丧失了固有的功能或者进化出了新的功能。

>Paralog: Paralogs are genes produced via gene duplication within a genome. Paralogs typically evolve new functions or else eventually become pseudogenes.

>**旁系／横向／并行同源（基因）**：指在一个特定的基因组中由于基因复制产生的同源。通常比较这些基因时，它们可能已经彼此具有了新的（不同）功能，也可能已经成为假基因了。

## 1.COG

**COG**是**C**luster of **O**rthologous **G**roups of proteins（蛋白相邻类的聚簇）的缩写。构成每个COG的蛋白都是被假定为来自于一个祖先蛋白，并且因此或者是orthologs或者是paralogs。

- **Orthologs**是指来自于不同物种的由垂直家系（物种形成）进化而来的蛋白，并且典型的保留与原始蛋白有相同的功能。
- **Paralogs**是那些在一定物种中的来源于基因复制的蛋白，可能会进化出新的与原来有关的功能。

**COG**是通过把所有完整测序的基因组的编码蛋白一个一个的互相比较确定的。在考虑来自一个给定基因组的蛋白时，这种比较将给出每个其他基因组的一个最相似的蛋白。这些基因的每一个都轮番的被考虑。如果在这些蛋白（或子集）之间一个相互的最佳匹配关系被发现，那么那些相互的最佳匹配将形成一个COG。这样，一个COG中的成员将与这个COG中的其他成员比起被比较的基因组中的其他蛋白更相像，尽管如果绝对相似性比较的。最佳匹配原则的使用，没有了人为选择的统计切除的限制，这就兼顾了进化慢和进化快的蛋白。然而，还有一个加的限制就是一个COG必须包含来自于3个种系发生上远的基因组的一个蛋白。

COG数据库主要包括细菌和真菌的蛋白信息；真核生物的一般称为KOG数据库。通过比对可以将某个蛋白序列注释到某一个COG中。根据直系同源序列构成的COG，如果序列近似，则可以推测该序列的功能。

### COG的分类

COG数据库可以根据特性分为3个方向，按照功能一共可以分为二十六类。

INFORMATION STORAGE AND PROCESSING

- [J] Translation, ribosomal structure and biogenesis
- [A] RNA processing and modification
- [K] Transcription
- [L] Replication, recombination and repair
- [B] Chromatin structure and dynamics

CELLULAR PROCESSES AND SIGNALING

- [D] Cell cycle control, cell division, chromosome partitioning
- [Y] Nuclear structure
- [V] Defense mechanisms
- [T] Signal transduction mechanisms
- [M] Cell wall/membrane/envelope biogenesis
- [N] Cell motility
- [Z] Cytoskeleton
- [W] Extracellular structures
- [U] Intracellular trafficking, secretion, and vesicular transport
- [O] Posttranslational modification, protein turnover, chaperones
- [X] Mobilome: prophages, transposons

METABOLISM

- [C] Energy production and conversion
- [G] Carbohydrate transport and metabolism
- [E] Amino acid transport and metabolism
- [F] Nucleotide transport and metabolism
- [H] Coenzyme transport and metabolism
- [I] Lipid transport and metabolism
- [P] Inorganic ion transport and metabolism
- [Q] Secondary metabolites biosynthesis, transport and catabolism

POORLY CHARACTERIZED

- [R] General function prediction only
- [S] Function unknown

最早从1997年开始推出COG，到2003年发布版本，2014年更新过一次，2020年由有一个新的版本。目前项目主页：https://www.ncbi.nlm.nih.gov/research/cog，可以在线检索获得COG信息以及携带该COG的细菌基因组Assembly信息，并链接到Genbank的蛋白质序列数据。

如果开发web应用，NCBI也可提供了相应的Web Services。如果要获得数据库信息，可以到 https://www.ncbi.nlm.nih.gov/research/cog-project/，访问NCBI Ftp下载。COG的软件可以到ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/COGsoft下载。

```bash
# 2003版 COGs
$ wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/myva
$ makeblastdb -in myva -dbtype prot -title cog -parse_seqids \
> -out PATH/TO/DB/cog -logfile PATH/TO/DB/cog.log

# 2014版 COGs
$ wget -r ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/* .
$ gunzip prot2003-2014.fa.gz
```

**2020版 COGs**

ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data

- cog-20.cog.tsv COG蛋白信息
- cog-20.def.tab COG注释信息
- cog-20.org.csv Genome Assembly
- cog-20.tax.csv 物种分类信息
- fun.tab COG分类信息

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

COG软件下载地址


## 参考资料

1. [NCBI的COG介绍](https://www.plob.org/article/1000.html)
2. [Expanded microbial genome coverage and improved protein family annotation in the COG database](https://www.ncbi.nlm.nih.gov/pubmed/25428365)
3. Fitch WM (1970) "Distinguishing homologous from analogous proteins". Syst Zool, 19:99-113.
4. Li WH, Yang J, Gu X (2005) "Expression divergence between duplicate genes". Trends Genet, 21:602-607.
5. http://blog.sciencenet.cn/blog-217859-280960.html
6. http://homepage.usask.ca/~ctl271/857/def_homolog.shtml
7. https://www.jianshu.com/p/f4060461c951
