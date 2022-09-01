# 新型冠状病毒基因组纳米孔测序常规流程

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! abstract "内容简介"
    本实验流程是参考 Sars-Cov-2 基因组纳米孔测序常见方案，根据我们实验室应对疫情和监测工作的实际流程，结合网络上的相关资源而建立的实验方案。本实验流程适合针对少量样本并希望获得较高通量数据的新冠基因组测序目标。

## 1. 实验方案简介

新冠病毒培养需要P3实验室，因此普通实验室只能针对样品核酸进行测序，如何把新冠病毒核酸富集到能够测序的水平是首先要解决的，虽然目前也有探针捕获等试剂盒或服务，但实验流程长，对于低浓度样本效果一般，不太适合疫情处置是快速溯源的需求。因此目前简便实用的方法是采用PCR扩增的方式。

本流程采用的方法是随机引物逆转录后多重PCR覆盖全基因组。这种方法也称为 Multiplex tilling PCR。扩增新冠病毒的引物方案较多，如ARTIC Network的实验方案采用短片段（400~500bp）进行扩增，有利于扩增获得样品中病毒序列较高的丰度。实验室可以根据病毒突变位点的变化选择引物方案并及时更新扩增引物版本。

采用Nanopore测序方法，对CT32以下的临床标本容易获得覆盖度较为完整的基因组数据。浓度越高越容易获得高质量的基因组序列。本方案适合多样本上机测序，如果样品数量较少，可以适当提高每个样本的投入量。

!!! info "测序质量"
    对于担心三代测序数据的质量低，是否会导致基因组不准确的问题。实际上我们对纳米孔测序结果突变位点进行过一代测序和二代测序的复核，结果均一致。虽然纳米孔测序目前的准确度约95%，但通过扩增子深度的叠加，错误碱基会被掩盖，只要具有一定的测序深度，一致性序列的准确度达到99.99%是没有问题的。因此利用artic标准分析流程，可以获得非常高质量的病毒一致性序列结果。因此我们认为测序和basecalling导致的错误几乎可以忽略，更应关注的是consensus序列分析过程中，软件的错误判别，特别是indel的判别，容易收到低覆盖度的其他序列污染或者多重感染等因素而误判。
    测序方法汇总可参考CDC的资料 https://github.com/CDCgov/SARS-CoV-2_Sequencing

## 2. 实验流程

### 2.1 cDNA获取

实验室目前针对送样都是采用自动化病毒核酸提取仪获得（我们使用的基本上均为国产品牌，如硕世、中元等），提取仪一般采用的都是磁珠法。如果使用手工提取试剂，如Qiagen等试剂盒，获得的产率可能会更好。不论那种方法，目前对于CT32以下的样本进行基因组测序获得完整基因组一般是没有问题的。对于CT32以上的样本，需要优化实验条件，才有可能获得较高覆盖度的基因组序列。对于CT35以上的样本，测序效果较差，一般不建议开展全基因组测序。

#### 2.1.1 试剂与耗材

- [NEB ]()
- [Invitrogen SuperScript™ IV First-Strand Synthesis System](https://www.thermofisher.com/order/catalog/product/18091050#/18091050)
- 0.2mL PCR管

!!! info "注意事项"
    SSIV的酶无RNaseH活性，获得长片段cDNA效果更好。SSIII，有中等活性的RNaseH活性，获得的cDNA进行短片段测序也可以成功，如果是长度更长一些的片段，SSIV效果明显更好一些。[不同酶的关键逆转录性能比较](https://upload-images.jianshu.io/upload_images/18923961-779adbe79350150a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
    在实际疫情处置过程中，对于时效性要求非常高，每一步都要尽可能节约时间，因此我们选择使用速度更快的NEB LunaScript预混试剂盒，其内含随机引物，逆转录实验只需10min即可完成。

#### 2.1.2 实验步骤

- 采用RNA试剂盒手工提取样本总RNA或使用提取仪自动提取。

>RNA通过Qiagen RNeasy Mini Kit 手工提取试剂盒或者机器提取的都可以，PCR实验结果CT无明显差异。有[文章](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6733840/)比较过Qiagen不同的RNA提取试剂盒的宏基因组测序reads产出效果，可以供大家参考。

- 取EP管，依次装入以下试剂，分装到0.2mL PCR中，加入RNA模板，PCR管用移液器吹吸混匀，离心将管壁残液甩下。如果需要增加使用的扩增引物数量，可以等体积增加混合成分。

**LunaScript逆转录**

| 试剂      | 1x        | 8x        |
| :-------- | :--------:| :--------:|
| LunaScript  | 2 µl | 16µl |
| Template RNA | 8 µl | 8µl * 8 |
| **Total** | **10 µl** | **10µl * 8** |

- 反应体系在PCR仪器中 25°C 2min。然后立即置于冰上至少1min。
- 反应体系在PCR仪器中 65°C 5min。然后立即置于冰上至少1min。

**SSIV逆转录**

| 试剂      | 1x        | 8x        |
| :-------- | :--------:| :--------:|
| 50µM random hexamers  | 1 µl | 8µl |
| 10mM dNTPs mix (10mM each) | 1 µl | 8µl |
| Template RNA | 11 µl | 11µl * 8 |
| **Total** | **13 µl** | **13µl * 8** |

!!! info "注意事项"
    1. 样本RNA做过SARS-CoV-2 qPCR测试，可以根据CT值初步判断RNA含量。如果RNA量在CT18～32之间，可以直接加入模板RNA。如果在12～15之间，100倍稀释；如果在15～18之间，10倍稀释。对于临床病例来说，大部分样本CT值在18～32之间。
    2. 使用SSIV random hexamers是为了最大程度获得全长cDNA，SSIV自带的random hexamers浓度是50ng/µl，6mers长度，如果加1µl，相当与摩尔数 25e-6µmol，因此根据表中量，需要加2µl。但根据试剂盒说明书中也是加1µl量，因此这里参考说明书。LunaScript MasterMix包含了6mer的随机引物，无需额外添加，如果要做特异性逆转录，则需要选择LunaScript其他试剂盒。
    3. 推荐加入样本前，在PCR仪上预热。

- 反应体系在PCR仪器中 65°C 5min。然后立即置于冰上至少1min。
- 按下表配置反应体系，添加到冰上的PCR管中。置于PCR仪器中，42°C 50min，70°C 10min，结束后5°C 保温。这一步要设置105°C热盖。

| 试剂 | 1x | 8x |
| :-------- | :--------: | :--------: |
| 5x SSIV Buffer | 4 µl | 32 µl |
| 100mM DTT| 1 µl | 8 µl |
| RNaseOUT RNase Inhibitor | 1 µl | 8 µl |
| SSIV Reverse Transcriptase | 1 µl | 8 µl |
| 前一步PCR管内液体 | 13 µl | 13 µl * 8 |
| **Total** | **20 µl** | **20 µl * 8** |

>以上配液操作在PCR工作站或生物安全柜内进行。使用前紫外处理30分钟并通风。

### 2.2 PCR扩增

目前实验室首选引物为Artic Network v4.1（针对奥密克戎优化，以下简称AN），Midnight（以下简称MD），VarSkip v1b（针对奥密克戎优化，以下简称VS）。由于国内大规模的疫情较少，一般处置的样本数量都是12个以下，而对于序列准确性要求高，因此常采用策略是每个样本同时使用2套及以上的引物进行扩增并独立测序，综合分析数据

目前使用的各组引物设计方案：

- artic network v4.1
- midnight
- varskip v
- varskip long

> 对于ligation法建立文库，本实验室一般采用AN v4.1和midnight引物组互相验证的方法，对于奥密克戎变异株具有较好的覆盖度。而varskip v1a在基因组26～27k处与AN v4.1

#### 2.2.1 试剂与耗材

- [Q5 Hot Start High-Fidelity DNA Polymerase](https://international.neb.com/products/m0493-q5-hot-start-high-fidelity-dna-polymerase)

>ARTIC Network流程里使用的是NEB Q5酶，试用Qiagen的多重PCR产品，效果也可以。

#### 2.2.2 实验步骤

这里以AN引物为例，其他引物参照配置：

- 将98组冻干引物 10000rpm 离心 5min。将每个引物配置成100uM浓度。
- 准备2个EP管，标记为P1和P2。将引物分成奇数组和偶数组，将配置的100uM浓度的引物，每管吸取5uL。取奇数组上下游引物到Pool1中，偶数组上下游引物到Pool2中。
- 震荡混匀，离心甩干后，各取50uL P1/2 引物，加入450uL不含核酸酶水，将Pool1和Pool2引物终浓度配置为10uM。

>以上配液操作在PCR工作站或生物安全柜内进行。使用前紫外处理30分钟并通风。

- 每个样本每组引物分别取2个PCR管（如同时做3组引物，则每个样本准备6个PCR管），使用NEB Q5 Hot Start DNA Polymerase试剂进行PCR扩增，配液参见下表。

根据引物区分配液，其中AN引物每管加3.6uL(或4uL)，VS引物加3uL，MD加1.1uL。

| NEB Q5 Hot Start DNA Polymerase | Pool1 | Pool2 |
| :-------- | :--------:| :--------:|
| 2X Q5 Master Mix | 12.5 µl | 12.5 µl |
| Primer Pool 1 or 2 (10µM) | 3.6 µl | 3.6 µl |
| Nuclease-free water |  6.4 µl | 6.4 µl |
| cDNA | 2.5uL | 2.5uL |
| **Total** | **25 µl** | **25 µl** |

>根据体系优化，每个引物的终浓度为0.015uM。

- NEB Q5 的PCR程序为：

| 序号 | 温度 | 时间 | 前往 |
| :-------- | :--------:| :--------:| :--------:|
| 1. | 98°C | 30s | |
| 2. | 98°C | 15s | |
| 3. | 63°C | 5min | 前往2, 循环25～35cycles |
| 4. | 4°C | Hold | |

>Ct18-21的样本用25循环, Ct>30的样本用35个循环。CT值介于21～30的根据分布自行设置。

>应用本方法到其他病原检测时，需要设计多重PCR的话，可以使用[这里](http://primal.zibraproject.org)的引物设计工具获得引物组。

### 2.3 Amplicon纯化回收

#### 2.3.1 试剂与耗材

- [Elution Buffer](https://www.qiagen.com/us/products/discovery-and-translational-research/lab-essentials/buffers-reagents/buffer-eb)
- [AMPure XP Beads](https://www.beckman.com/reagents/genomic/cleanup-and-size-selection/pcr)

#### 2.3.2 实验步骤

- 将每个样本的Pool1和Pool2扩增产物混合，因为测序文库与上样量有关，本实验流程设计的是针对8个以上样本的上样量。建议检测样本带质控DNA做为阴性对照，以便检查是否有污染序列。
- 将AMPure XP Beads 提前30min取出置于常温，使用前在震荡仪上混合均匀。
- 取50uL AMPure XP Beads 到 50uL 混合pooling Amplicon中。震荡混匀（增加震荡时间能极大提高回收率），离心甩干。

> 1x 磁珠吸附250～1000bp效率较好。具体[参见](https://research.fhcrc.org/content/dam/stripe/hahn/methods/mol_biol/SPRIselect%20User%20Guide.pdf)。对于1200bp的MD引物扩增产物，也可以采用0.6x比例的磁珠进行回收。

- 室温放置5min。
- 将EP管放置于磁力架上2min，直至液体澄清。
- 吸弃液体，用200uL新鲜配置的70～80%乙醇洗2次。
- 用10uL移液器吸弃残液，开盖1min，让乙醇挥发。
- 从磁力架上取下，用30uL Elution Buffer重悬磁珠，用枪头或手指轻轻混匀，静置2min。
- 置于磁力架上，液体澄清后吸取30uL液体到一个新的EP管中。取1uL进行Qubit定量。

>对于室温放置，也可置于37°C金属浴。延长温度及37°C处置，更适合长片段DNA的回收。

### 2.4 文库构建

本方案采用ligation法进行文库构建。

#### 2.4.1 试剂与耗材

- Nanopore EXP-NBD104/114
- Nanopore SQK-LSK109
- Nanopore 测序芯片制备试剂盒 EXP-FLP002
- Eppendorf DNA LoBind Tubes
- [NEBNext® Ultra™ II DNA Library Prep Kit for Illumina E7645](https://international.neb.com/products/e7645-nebnext-ultra-ii-dna-library-prep-kit-for-illumina)

#### 2.4.2 实验步骤

>注意此流程是根据6~24个样本混样的量规划的。如果只做单个样本，可以不加barcodes，但需要调整input DNA量，建议加入DNA达40ng，保证上机时能约有20ng DNA，达到最佳的芯片测序分子数。

- 将定量的DNA稀释成1ng/uL，该浓度是针对700bp长度的Amplicon，如果片段长400bp，建议浓度为2ng/uL，这样分子浓度约为100～200fmol。
- 取0.2mLPCR管，根据下表进行配液，进行末端修复。

| 试剂 | 量 |
| :----- | :--: |
| DNA amplicons | 5 µl |
| Nuclease-free water | 7.5 µl |
| Ultra II End Prep Reaction Buffer | 1.75 µl |
| Ultra II End Prep Enzyme Mix | 0.75 µl |
| **总量** | **15 µl** |

- 室温放置10min。置于PCR仪上，65°C 5min，立即置于冰上至少1min
- 在PCR管中加入以下试剂：

| 试剂 | 量 |
| :----- | :--: |
| 连接NBXX的混合液 | 15 µl |
| NBXX barcode | 2.5 µl |
| Ultra II Ligation Master Mix | 17.5 µl |
| Ligation Enhancer | 0.5 µl |
| **总量** | **35.5 µl** |

>NBXX barcodes 为 EXP-NBD104(01~12)和EXP-NBD114(13～24) 试剂中的 barcodes

- 室温放置15min。置于PCR仪上，70°C 10min，立即置于冰上至少1min

>70°C 10min 目的是抑制DNA Ligase活性，避免接头形成二聚体。

- 将连接barcodes后的所有样本的混合到一个EP管中。取1uL进行定量
- 取LoBind Tubes EP管，按照下表进行配液，室温放置15min。

| 试剂 | 量 |
| :----- | :--: |
| Barcoded amplicon pools | 30 µl |
| NEBNext Quick Ligation Reaction Buffer (5X) | 10 µl |
| AMII adapter mix | 5 µl |
| Quick T4 DNA Ligase | 5 µl |
| **总量** | **50 µl** |

>input DNA的量根据barcods数量决定，数量在40ng～160ng(8～24 barcods)。

- 加入50uL AMPure XP Beads，轻轻混匀。室温放置5min
- 置于磁力架上2min，液体澄清后吸弃液体。
- 磁力架上取下EP管，加入200uL SFB，重悬磁珠。
- EP管放到磁力架上静置2min，液体澄清后吸弃。用SFB再洗1次。
- 加入15uL EB重悬磁珠，室温放置2min。
- 将EP管置于磁力架上。吸取澄清液体到一个新的EP管中，为测序文库。

>用于上机测序的文库DNA量为20ng

- 将30uL FLT 加入到1管FB中，震荡混匀。
- 打开测序芯片的priming port，用1mL移液器垂直顶住priming port，慢慢旋转移液器量程，黄色纳米孔保护液，直到没有气泡。
- 吸取800uL FLT/FB混合液，从Priming port中加入，可以不用加完，避免气泡加入芯片中。等待5min。

![芯片Priming port和SpotON port Cover位置介绍](https://upload-images.jianshu.io/upload_images/18923961-6970fac21978e41e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 打开Spot on port，滴加200uL FLT/FB混合液。直至Spoton port完全吸收。
- 取一个EP管，按照下表配液：

| 试剂 | 量 |
| :----- | :--: |
| SQB | 37.5 µl |
| LB | 25.5 µl |
| Final library | 12 µl |
| **总量** | **75 µl** |

>LB使用前再充分混匀

- 滴加75uL文库到Spoton Port中，液体完全吸收后关闭Spoton Port和Priming Port，将芯片装入测序仪，打开MinKnow测序软件进行测序。

>对于没用用完nanopores的芯片，可以启用清洗程序以后继续使用：

- 取 1管 Wash Solution B置于室温，震荡混匀后放在冰上待用。
- 在 1个 EP管中，加入20uL Wash Solution A和380uL Wash Solution B。吹吸混匀（不要震荡），置于冰上。
- 暂停测序，但不要取出芯片。
- 确定Priming Port和SpotON Port关闭。用1000uL移液器，从Waste Port1中吸出所有液体，确定芯片纳米孔处没有残液。
- 打开Priming Port，先确认消除气泡。然后加入400Wash Solution洗液，注意避免加入气泡。
- 关闭Priming Port，等30min后，用1000uL移液器从Waste Port1中吸弃所有液体。
- 如果还要测其他样本，可参考上样流程操作。如果暂时不用，则按照下面流程进行保存芯片。
- 将一管Storage Buffer(S)取出置于室温，溶解后上下颠倒混匀。
- 打开用完的芯片Priming Port，1000uL移液器设置成200uL量程，然后枪头Priming Port顶住孔，旋转移液器量程20uL，以确定排除气泡。
- 吸取500uL Storage Buffer(S)，从Priming Port中缓慢加入。关闭Priming Port。
- 关闭Priming Port和SpotON Port，用1000uL移液器从Waste Port1中吸取所有液体。
- 芯片置于冰箱冷藏。

>下一次使用的芯片，最好设置成不同的电压以获得更好的结果，[参见](https://community.nanoporetech.com/protocols/experiment-companion-minknow/v/mke_1013_v1_revar_11apr2016/adjusting-the-starting-pot)

## 3. 数据分析

### 3.1 软件安装

使用 conda 构建独立的运行环境。

```bash
$ git clone --recursive https://github.com/artic-network/artic-ncov2019.git
$ conda env create -n ncov -f artic-ncov2019/environment.yml
$ conda activate ncov
(ncov)$ conda list
```

如果使用docker images来运行guppy

```bash
$ docker pull genomicpariscentre/guppy
# 如果服务器有GPU支持，可以选择guppy-gpu
$ docker pull genomicpariscentre/guppy-gpu

# 交互模式运行container
$ docker run --name my_exp -it -v local_dir_of_fast5:/media/ genomicpariscentre/guppy /bin/bash

# 对fast5数据进行basecalling
$ guppy_basecaller -i $HOME/data/fast5 -s output_folder \
> -c dna_r9.4.1_450bps_hac.cfg --barcode_kits EXP-NBD104 \
> --num_callers 40 --trim_barcodes
```

### 3.2 分析流程

> 与illumina原始数据是图像文件不同，Nanopore的原始数据以HDF5格式保存。HDF是"Hierarchical Data Format"首字母缩写，HDF5是一种纯文本，类似json的数据格式记录电信号。查看HDF5原始数据可以用 [hdf5_tools](https://support.hdfgroup.org/products/hdf5_tools/)工具。数据文件后缀一般用.fast5表示。

!!! info "背景知识"
    Nanopore 测序数据期初由于其错误率高而“闻名”，特别是对二聚体的区分。电信号数据随着处理软件和算法的不断改进，准确率得到不断的提升。HDF5数据进行basecalling的软件目前很多，采用的算法也各异。官方的如albacore，guppy等，第三方工具也有很多。但随着官方guppy的出现和RNN算法的引入，目前测序数据基本上只需要采用guppy进行basecalling即可。对于模型算法推荐super accuracy model，但这个算法需要更强大的GPU支持，即使在Gridion的GV100上，大概也只能同时跑满3张flow cells。对于例如bonito等软件，虽然可能准确度能更进一步提升，但消耗的计算资源也更高。

Gridion机器上的MinKnow软件目前默认使用guppy进行basecalling，目前版本(v5.0.17)支持demultiplex。默认生成的fastq可直接用于分析。如果需要采用不同的软件提升质量，可以将fast5_pass中的fast5数据复制到其他服务器，再进行basecalling。

```bash
# 在服务器上 rsync 同步数据到本地
(ncov)$ rsync -avz username@gridion:/path/to/fast5 .
# 或者在gridion上同步数据到服务器端
$ rsync -avz /path/to/fast5 username@server:~/data/fast5
#这里设置服务器端 fast5 数据目录存放于：$HOME/data/fast5；采用的是R9.4的测序芯片，各个不同测序芯片basecalling设置参数可以参见[这里](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html)

# guppy sup 模式basecalling
(ncov)$ guppy_basecaller -c dna_r9.4.1_450bps_sup.cfg \
> -i $HOME/data/fast5 -s fastq_output -x 'auto'
```

### 3.3 手动分析流程

```bash
# 下载参考基因组序列
$ efetch -id NC_045512 -db nuccore -format fasta > NC_045512.fasta
# 生成样本的 fastq
$ for i in 'fastq_pass/*.fastq'; do gunzip -c $i | gzip >> sample_1.fq.gz; done
# 去除可能的接头序列
$ porechop -t 4 -i sample_1.fq.gz -o sample_1.clean.fq.gz
# 去除低质量序列
$ gunzip -c sample_1.clean.fq.gz | NanoFilt -q 10 -l 400 --maxlength 700 | gzip > sample_1.clean.filt.fq.gz
# 测序数据与参考基因组进行比对，生成bam比对结果
$ bwa index NC_045512.fasta
$ bwa mem -R '@RG\t\ID:sample1\tSM:sample1' -t 4 -x ont2d NC_045512.fasta sample_1.clean.filt.fq.gz | \
> samtools view -bS | \
> samtools sort -o sample_1.clean.filt.sorted.bam
$ samtools index sample_1.clean.filt.sorted.bam
$ samtools faidx NC_045512.fasta
# polishing and vcf calling
$ medaka consensus

$ bcftools mpileup -f --threads 4 MN908047.3.fasta sample_1.clean.filt.sorted.bam | \
> bcftools call --ploidy 1 -mv -Oz > sample_1.vcf.gz
$ tabix sample_1.vcf.gz
# 获得 consensus 序列
$ bcftools consensus -f MN908947.3.fasta -H 1 \
> sample_1.vcf.gz -o sample_1.consensus.fasta

# 或者使用minimap工具进行比对
$ minimap2
```

## 4. 参考资料

- https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w
- https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
- https://nanoporetech.com/resource-centre/nanopore-sequencing-sars-cov-2-genome-introduction-protocol
