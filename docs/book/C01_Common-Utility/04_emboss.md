# EMBOSS

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

![emboss](http://emboss.sourceforge.net/images/emboss.jpg)

!!! Abstract "内容简介"
    EMBOSS 是欧洲分子生物学组织开发的 Unix/Linux 下的生物学分析工具。EMBOSS 包含工具众多，这里只介绍与微生物基因组分析可能会用到的一些工具，所有的软件和其文档参考官方文档。

由于 EMBOSS 的命令行方式比较传统，类似phylip的交互模式，和当前许多高通量测序工具的参数设置等不一致，加上许多新的工具代替，因此使用的机会不是很多。这里介绍部分与微生物基因组分析可能会用到的一些工具，所有的软件和其文档参考[官方文档](http://emboss.sourceforge.net/docs/)。

## 1. 安装 EMBOSS

```bash
# ubuntu 包含 emboss 发行版
$ sudo apt install emboss

# 通过 conda 安装
$ conda create -n emboss emboss
$ conda activate emboss
```

EMBOSS 工具集包含众多基于命令行的工具，可以集成到分析工具流中。EMBOSS 中的许多可以直接访问远程数据库，但这需要默认配置。如果是通过 conda 安装，则要将所在conda虚拟环境路径中的`share/EMBOSS/emboss.default.template`复制为`emboss.default`。如果是通过系统级安装，要在当前用户环境中配置数据库，则可以将`emboss.default.template`复制为`~/.embossrc`。

```bash
# 建立配置文件
(emboss)$ cd $CONDA_PREFIX/EMBOSS
(emboss)$ cp emboss.default.template emboss.default
(emboss)$ vim emboss.default
# 将需要激活的数据库注释符号取出，比如embl

# 查看可以使用的数据库
(emboss)$ showdb
```

## 2. 使用

### 查看各个程序

**wossname**: EMBOSS 包含众多应用程序，为方便查询，有一个专门的程序`wossname`可以查看所有程序及功能。

```bash
# 显示所有程序及功能
(emboss)$ wossname
# 显示某一个类别的程序
# 显示alignment比对相关的应用程序
(emboss)$ wossname alignment

# 按功能分组显示应用
(emboss)$ wossname -search "" | less

# 按字母顺序显示应用
(emboss)$ wossname -search "" -alphabetic | less
```

### 数据库相关

**showdb**: 显示可用的数据库
**dbtell**: 显示数据库相关信息

```bash
# 显示当前可用数据库
(emboss)$ showdb -full

# 显示embl数据库信息
(emboss)$ dbtell embl
```


### 序列处理

**[seqret](http://structure.usc.edu/emboss/seqret.html)**: 生成序列或格式化序列

```bash
# 读入序列，以60字符为长度将其格式化为规范的fasta格式数据
(emboss)$ seqret 1.fas 2.fas

# 序列个是转换 fasta 格式转换为其他格式
(emboss)$ seqret fasta::1.fas gcg::1.gcg

# 直接从服务器读取序列
(emboss)$ seqret embl:x52524

# 通过交互方式进行操作
(emboss)$ seqret
Read and write (return) sequences
Input (gapped) sequence(s):         <--- 输入序列文件名称
output sequence(s) [...]:           <--- 输入序列保存的文件名称
```

**[showfeat](http://structure.usc.edu/emboss/showfeat.html)**: 显示序列信息

```bash
# 显示组装子各个 contigs 长度
(emboss)$ showfeat assembly.fasta result.showfeat
(emboss)$ less result.showfeat
```

**[transeq](http://structure.usc.edu/emboss/transeq.html)**: DNA/RNA序列翻译成氨基酸序列，文件以 `.pep` 后缀保存  
**[backtranseq](http://structure.usc.edu/emboss/backtranseq.html)**: 氨基酸序列转换成DNA序列。

```bash
# fasta个是的DNA序列转换成氨基酸序列
(emboss)$ transeq input.fasta output.pep
# 将 input.pep 氨基酸序列转换成碱基格式保存成 output.fas DNA序列
(emboss)$ backtranseq input.pep output.fas
```

### 引物设计

**eprimer3**:
**primersearch**:

### 双序列比对

双序列(Pairwise Alignment)比对的软件

**needle**: local alignment 序列比对
**water**: global alignment 序列比对，基于 water-smith 算法。

```bash
# 对2条序列 seq1.fas 和 seq2.fas 进行全局比对，生成 seq1v2 alignment
(emboss)$ water seq1.fas seq2.fas seq1v2.sw
(emboss)$ needle seq1.fas seq2.fas seq1v2.nl
```

**needleall**

### 多序列比对

**emma**: 基于 clustalW 的多序列比对程序
**edialign**: Local 多序列比对

```bash
(emboss)$ emma msa1.fas msa1.phy
(emboss)$ edialign
```

### 数据可视化

**density**: 绘制核酸密度图。

```bash
# 生成序列密度图
(emboss)$ density -seqall input.fasta -display D -gragh ps
```

## 工具列表

具体信息参见：[工具列表](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/)

| 程序名称 | 用途 |
| -------- | ---- |
| aaindexextract | 从AAINDEX中提取氨基酸属性 |
| abiview | 显示 ABI 测序数据峰值 |
| acdc | 测试 ACD 文件 |
| acdpretty | 修正 ACD 格式文件 |
| acdtable | 从 ACD 文件应用中生成 HTML 格式的表格 |
| acdtrace | Trace processing of an application ACD file (for testing) |
| acdvalid | Validate an application ACD file |
| aligncopy | Read and write alignments |
| aligncopypair | Read and write pairs from alignments |
| antigenic | Find antigenic sites in proteins |
| assemblyget | Get assembly of sequence reads |
| backtranambig | Back-translate a protein sequence to ambiguous nucleotide sequence |
| backtranseq | Back-translate a protein sequence to a nucleotide sequence |
| banana | Plot bending and curvature data for B-DNA |
| biosed | Replace or delete sequence sections |
| btwisted | Calculate the twisting in a B-DNA sequence |
| cachedas | Generate server cache file for DAS servers or for the DAS registry |
| cachedbfetch | Generate server cache file for Dbfetch/WSDbfetch data sources |
| cacheebeyesearch | Generate server cache file for EB-eye search domains |
| cacheensembl | Generate server cache file for an Ensembl server |
| cai | 计算 condon adaptation index |
| chaos | Draw a chaos game representation plot for a nucleotide sequence |
| charge | Draw a protein charge plot |
| checktrans | Report STOP codons and ORF statistics of a protein |
| chips | Calculate Nc codon usage statistic |
| cirdna | Draw circular map of DNA constructs |
| codcmp | Codon usage table comparison |
| codcopy | Copy and reformat a codon usage table |
| coderet | Extract CDS, mRNA and translations from feature tables |
| compseq | Calculate the composition of unique words in sequences |
| cons | Create a consensus sequence from a multiple alignment |
| consambig | Create an ambiguous consensus sequence from a multiple alignment |
| cpgplot | Identify and plot CpG islands in nucleotide sequence(s) |
| cpgreport | Identify and report CpG-rich regions in nucleotide sequence(s) |
| cusp | Create a codon usage table from nucleotide sequence(s) |
| cutgextract | Extract codon usage tables from CUTG database |
| cutseq | Remove a section from a sequence |
| dan | Calculate nucleic acid melting temperature |
| dbiblast | Index a BLAST database |
| dbifasta | Index a fasta file database |
| dbiflat | Index a flat file database |
| dbigcg | Index a GCG formatted database |
| dbtell | Display information about a public database |
| dbxcompress | Compress an uncompressed dbx index |
| dbxedam | Index the EDAM ontology using b+tree indices |
| dbxfasta | Index a fasta file database using b+tree indices |
| dbxflat | Index a flat file database using b+tree indices |
| dbxgcg | Index a GCG formatted database using b+tree indices |
| dbxobo | Index an obo ontology using b+tree indices |
| dbxreport | Validate index and report internals for dbx databases |
| dbxresource | Index a data resource catalogue using b+tree indices |
| dbxstat | Dump statistics for dbx databases |
| dbxtax | Index NCBI taxonomy using b+tree indices |
| dbxuncompress | Uncompress a compressed dbx index |
| degapseq | Remove non-alphabetic (e.g. gap) characters from sequences |
| density | Draw a nucleic acid density plot |
| descseq | Alter the name or description of a sequence |
| diffseq | Compare and report features of two similar sequences |
| distmat | Create a distance matrix from a multiple sequence alignment |
| dotmatcher | Draw a threshold dotplot of two sequences |
| dotpath | Draw a non-overlapping wordmatch dotplot of two sequences |
| dottup | Display a wordmatch dotplot of two sequences |
| dreg | Regular expression search of nucleotide sequence(s) |
| drfinddata | Find public databases by data type |
| drfindformat | Find public databases by format |
| drfindid | Find public databases by identifier |
| drfindresource | Find public databases by resource |
| drget | Get data resource entries |
| drtext | Get data resource entries complete text |
| edamdef | Find EDAM ontology terms by definition |
| edamhasinput | Find EDAM ontology terms by has_input relation |
| edamhasoutput | Find EDAM ontology terms by has_output relation |
| edamisformat | Find EDAM ontology terms by is_format_of relation |
| edamisid | Find EDAM ontology terms by is_identifier_of relation |
| edamname | Find EDAM ontology terms by name |
| edialign | Local multiple alignment of sequences |
| einverted | Find inverted repeats in nucleotide sequences |
| embossdata | Find and retrieve EMBOSS data files |
| embossupdate | Checks for more recent updates to EMBOSS |
| embossversion | Report the current EMBOSS version number |
| emma | Multiple sequence alignment (ClustalW wrapper) |
| emowse | Search protein sequences by digest fragment molecular weight |
| entret | Retrieve sequence entries from flatfile databases and files |
| epestfind | Find PEST motifs as potential proteolytic cleavage sites |
| eprimer3 | Pick PCR primers and hybridization oligos |
| eprimer32 | Pick PCR primers and hybridization oligos |
| equicktandem | Find tandem repeats in nucleotide sequences |
| est2genome | Align EST sequences to genomic DNA sequence |
| etandem | Find tandem repeats in a nucleotide sequence |
| extractalign | Extract regions from a sequence alignment |
| extractfeat | Extract features from sequence(s) |
| extractseq | Extract regions from a sequence |
| featcopy | Read and write a feature table |
| featmerge | Merge two overlapping feature tables |
| featreport | Read and write a feature table |
| feattext | Return a feature table original text |
| findkm | Calculate and plot enzyme reaction data |
| freak | Generate residue/base frequency table or plot |
| fuzznuc | Search for patterns in nucleotide sequences |
| fuzzpro | Search for patterns in protein sequences |
| fuzztran | Search for patterns in protein sequences (translated) |
| garnier | Predict protein secondary structure using GOR method |
| geecee | Calculate fractional GC content of nucleic acid sequences |
| getorf | Find and extract open reading frames (ORFs) |
| godef | Find GO ontology terms by definition |
| goname | Find GO ontology terms by name |
| helixturnhelix | Identify nucleic acid-binding motifs in protein sequences |
| hmoment | Calculate and plot hydrophobic moment for protein sequence(s) |
| iep | Calculate the isoelectric point of proteins |
| infoalign | Display basic information about a multiple sequence alignment |
| infoassembly | Display information about assemblies |
| infobase | Return information on a given nucleotide base |
| inforesidue | Return information on a given amino acid residue |
| infoseq | Display basic information about sequences |
| isochore | Plot isochores in DNA sequences |
| jaspextract | Extract data from JASPAR |
| jaspscan | Scan DNA sequences for transcription factors |
| lindna | Draw linear maps of DNA constructs |
| listor | Write a list file of the logical OR of two sets of sequences |
| makenucseq | Create random nucleotide sequences |
| makeprotseq | Create random protein sequences |
| marscan | Find matrix/scaffold recognition (MRS) signatures in DNA sequences |
| maskambignuc | Mask all ambiguity characters in nucleotide sequences with N |
| maskambigprot | Mask all ambiguity characters in protein sequences with X |
| maskfeat | Write a sequence with masked features |
| maskseq | Write a sequence with masked regions |
| matcher | Waterman-Eggert local alignment of two sequences |
| megamerger | Merge two large overlapping DNA sequences |
| merger | Merge two overlapping sequences |
| msbar | Mutate a sequence |
| mwcontam | Find weights common to multiple molecular weights files |
| mwfilter | Filter noisy data from molecular weights file |
| needle | Needleman-Wunsch 全局比对 |
| needleall | 两两双序列比对 |
| newcpgreport | Identify CpG islands in nucleotide sequence(s) |
| newcpgseek | Identify and report CpG-rich regions in nucleotide sequence(s) |
| newseq | Create a sequence file from a typed-in sequence |
| nohtml | Remove mark-up (e.g. HTML tags) from an ASCII text file |
| noreturn | Remove carriage return from ASCII files |
| nospace | Remove whitespace from an ASCII text file |
| notab | Replace tabs with spaces in an ASCII text file |
| notseq | Write to file a subset of an input stream of sequences |
| nthseq | Write to file a single sequence from an input stream of sequences |
| nthseqset | Read and write (return) one set of sequences from many |
| octanol | Draw a White-Wimley protein hydropathy plot |
| oddcomp | Identify proteins with specified sequence word composition |
| ontocount | Count ontology term(s) |
| ontoget | Get ontology term(s) |
| ontogetcommon | Get common ancestor for terms |
| ontogetdown | Get ontology term(s) by parent id |
| ontogetobsolete | Get ontology ontology terms |
| ontogetroot | Get ontology root terms by child identifier |
| ontogetsibs | Get ontology term(s) by id with common parent |
| ontogetup | Get ontology term(s) by id of child |
| ontoisobsolete | Report whether an ontology term id is obsolete |
| ontotext | Get ontology term(s) original full text |
| palindrome | Find inverted repeats in nucleotide sequence(s) |
| pasteseq | Insert one sequence into another |
| patmatdb | Search protein sequences with a sequence motif |
| patmatmotifs | Scan a protein sequence with motifs from the PROSITE database |
| pepcoil | Predict coiled coil regions in protein sequences |
| pepdigest | Report on protein proteolytic enzyme or reagent cleavage sites |
| pepinfo | Plot amino acid properties of a protein sequence in parallel |
| pepnet | Draw a helical net for a protein sequence |
| pepstats | Calculate statistics of protein properties |
| pepwheel | Draw a helical wheel diagram for a protein sequence |
| pepwindow | Draw a hydropathy plot for a protein sequence |
| pepwindowall | Draw Kyte-Doolittle hydropathy plot for a protein alignment |
| plotcon | Plot conservation of a sequence alignment |
| plotorf | Plot potential open reading frames in a nucleotide sequence |
| polydot | Draw dotplots for all-against-all comparison of a sequence set |
| preg | Regular expression search of protein sequence(s) |
| prettyplot | Draw a sequence alignment with pretty formatting |
| prettyseq | Write a nucleotide sequence and its translation to file |
| primersearch | Search DNA sequences for matches with primer pairs |
| printsextract | Extract data from PRINTS database for use by pscan |
| profit | Scan one or more sequences with a simple frequency matrix |
| prophecy | Create frequency matrix or profile from a multiple alignment |
| prophet | Scan one or more sequences with a Gribskov or Henikoff profile |
| prosextract | Process the PROSITE motif database for use by patmatmotifs |
| pscan | Scan protein sequence(s) with fingerprints from the PRINTS database |
| psiphi | Calculates phi and psi torsion angles from protein coordinates |
| rebaseextract | Process the REBASE database for use by restriction enzyme applications |
| recoder | Find restriction sites to remove (mutate) with no translation change |
| redata | Retrieve information from REBASE restriction enzyme database |
| refseqget | Get reference sequence |
| remap | Display restriction enzyme binding sites in a nucleotide sequence |
| restover | Find restriction enzymes producing a specific overhang |
| restrict | Report restriction enzyme cleavage sites in a nucleotide sequence |
| revseq | Reverse and complement a nucleotide sequence |
| seealso | Find programs with similar function to a specified program |
| seqcount | Read and count sequences |
| seqmatchall | All-against-all word comparison of a sequence set |
| seqret | Read and write (return) sequences |
| seqretsetall | Read and write (return) many sets of sequences |
| seqretsplit | Read sequences and write them to individual files |
| seqxref | Retrieve all database cross-references for a sequence entry |
| seqxrefget | Retrieve all cross-referenced data for a sequence entry |
| servertell | Display information about a public server |
| showalign | Display a multiple sequence alignment in pretty format |
| [showdb](http://structure.usc.edu/emboss/showdb.html) | 显示 EMBOSS 工具支持的数据库，一些工具可以直接操作远程数据库 |
| showfeat | 显示序列的属性，如长度等信息 |
| showorf | Display a nucleotide sequence and translation in pretty format |
| showpep | Display protein sequences with features in pretty format |
| showseq | Display sequences with features in pretty format |
| showserver | Display information on configured servers |
| shuffleseq | Shuffle a set of sequences maintaining composition |
| sigcleave | Report on signal cleavage sites in a protein sequence |
| silent | Find restriction sites to insert (mutate) with no translation change |
| sirna | Find siRNA duplexes in mRNA |
| sixpack | Display a DNA sequence with 6-frame translation and ORFs |
| sizeseq | Sort sequences by size |
| skipredundant | Remove redundant sequences from an input set |
| skipseq | Read and write (return) sequences, skipping first few |
| splitsource | Split sequence(s) into original source sequences |
| splitter | Split sequence(s) into smaller sequences |
| stretcher | Needleman-Wunsch rapid global alignment of two sequences |
| stssearch | Search a DNA database for matches with a set of STS primers |
| supermatcher | Calculate approximate local pair-wise alignments of larger sequences |
| syco | Draw synonymous codon usage statistic plot for a nucleotide sequence |
| taxget | Get taxon(s) |
| taxgetdown | Get descendants of taxon(s) |
| taxgetrank | Get parents of taxon(s) |
| taxgetspecies | Get all species under taxon(s) |
| taxgetup | Get parents of taxon(s) |
| tcode | Identify protein-coding regions using Fickett TESTCODE statistic |
| textget | Get text data entries |
| textsearch | Search the textual description of sequence(s) |
| tfextract | Process TRANSFAC transcription factor database for use by tfscan |
| tfm | Display full documentation for an application |
| tfscan | Identify transcription factor binding sites in DNA sequences |
| tmap | Predict and plot transmembrane segments in protein sequences |
| tranalign | Generate an alignment of nucleic coding regions from aligned proteins |
| transeq | Translate nucleic acid sequences |
| trimest | Remove poly-A tails from nucleotide sequences |
| trimseq | Remove unwanted characters from start and end of sequence(s) |
| trimspace | Remove extra whitespace from an ASCII text file |
| twofeat | Find neighbouring pairs of features in sequence(s) |
| union | Concatenate multiple sequences into a single sequence |
| urlget | Get URLs of data resources |
| variationget | Get sequence variations |
| vectorstrip | Remove vectors from the ends of nucleotide sequence(s) |
| water | Smith-Waterman local alignment of sequences |
| whichdb | Search all sequence databases for an entry and retrieve it |
| wobble | Plot third base position variability in a nucleotide sequence |
| wordcount | Count and extract unique words in molecular sequence(s) |
| wordfinder | Match large sequences against one or more other sequences |
| wordmatch | Find regions of identity (exact matches) of two sequences |
| wossdata | Find programs by EDAM data |
| wossinput | Find programs by EDAM input data |
| wossname | Find programs by keywords in their short description |
| wossoperation | Find programs by EDAM operation |
| wossoutput | Find programs by EDAM output data |
| wossparam | Find programs by EDAM parameter |
| wosstopic | Find programs by EDAM topic |
| yank | Add a sequence reference (a full USA) to a list file |
