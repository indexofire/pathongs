# 新型冠状病毒基因组测序

---

本节围绕目前新型冠状病毒基因组测序，根据我们自己实验室的工作流程结合网络上的资源，给需要开展病毒基因组测序的实验室提供帮助。

## 实验方法

目前常见的新型冠状病毒基因组测序方法较多，这里主要针对ARTIC Network的实验方案进行介绍。该方案采用分组的多重PCR扩增500bp左右的小片段覆盖全基因组，采用Nanopore测序方法进行。对CT30以下的临床标本提取RNA可以获得较好的基因组数据。CT30以上的样本，可能会出现部分区域gap的情况。我们采用该方法在illumina miseq平台也获得了较好的结果。实验方法一并在这里介绍。

> 测序方法汇总可参考CDC的资料 https://github.com/CDCgov/SARS-CoV-2_Sequencing

---

> 本文档采用 ARTIC Network 流程进行实验，可用于illumina平台和nanopore平台对SARS-Cov-2病毒进行基因组测序。

@indexofire <indexofire@gmail.com>

# 1. cDNA获取

## 1.1 试剂与耗材

- RNA提取试剂
- [Invitrogen SuperScript™ IV First-Strand Synthesis System](https://www.thermofisher.com/order/catalog/product/18091050#/18091050)
- 0.2mL PCR管

>SSIV的酶无RNaseH活性，获得长片段(～12kb)的全长cDNA效果更好。我们使用SSIII，有中等活性的RNaseH活性，获得的cDNA进行测序也可以成功，在目标reads丰度上SSIV似乎更好一些。

![不同酶的关键逆转录性能比较](https://upload-images.jianshu.io/upload_images/18923961-779adbe79350150a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 1.2 实验流程

- 采用RNA试剂盒手工提取样本总RNA，也可以用机器自动提取。

>RNA通过Qiagen RNeasy Mini Kit 手工提取试剂盒或者机器提取的都可以，PCR实验结果CT无明显差异。有[文章](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6733840/)比较过Qiagen不同的RNA提取试剂盒的宏基因组测序reads产出效果，可以供大家参考。

- 取EP管，依次装入以下试剂，分装到0.2mL PCR中，加入RNA模板，PCR管用移液器吹吸混匀，离心将管壁残液甩下。

| 试剂 | 1x | 8x |
| :-------- | :--------:| :--------:|
| 50µM random hexamers | 1 µl | 8µl |
| 10mM dNTPs mix (10mM each) | 1 µl | 8µl |
| Template RNA | 11 µl | 11µl * 8 |
| **Total** | **13 µl** | **13µl * 8** |

>注意事项：样本RNA做过SARS-CoV-2 qPCR测试，可以根据CT值初步判断RNA含量。如果RNA量在CT18～35之间，可以直接加入模板RNA。如果在12～15之间，100倍稀释；如果在15～18之间，10倍稀释。对于临床病例来说，大部分样本CT值在18～35之间。

>使用random hexamers是为了最大程度获得全长cDNA，SSIV自带的random hexamers浓度是50ng/µl，6mers长度，如果加1µl，相当与摩尔数 25e-6µmol，因此根据表中量，需要加2µl。但根据试剂盒说明书中也是加1µl量，因此这里参考说明书。
>推荐加入样本前，在PCR仪上预热。

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

# 2. PCR扩增

## 2.1 试剂与耗材

- DEPC H2O
- [QIAGEN Multiplex PCR Kit](https://www.qiagen.com/us/products/discovery-and-translational-research/pcr-qpcr/pcr-enzymes-and-kits/end-point-pcr/qiagen-multiplex-pcr-kit/)
- [Q5 Hot Start High-Fidelity DNA Polymerase](https://international.neb.com/products/m0493-q5-hot-start-high-fidelity-dna-polymerase)

>ARTIC Network流程里使用的是NEB Q5酶，试用Qiagen的多重PCR产品，效果也可以。

## 2.2 实验流程

- 将98组冻干引物 8000rpm 离心 10min。将每个引物配置成100uM浓度。
- 准备2个EP管，标记为P1和P2。从98个100uM的引物，每管吸取5uL分别到Pool1和2中，最终每管含100uM浓度的引物共980uL。
- 震荡混匀，离心甩干后，各取20uL P1/2 引物，加入180uL水。使得Pool1和Pool2引物终浓度为10uM。

>以上配液操作在PCR工作站或生物安全柜内进行。使用前紫外处理30分钟并通风。

- 每个样本取2个PCR管，使用NEB Q5 Hot Start DNA Polymerase试剂或者Qiagen Multiplex PCR试剂进行PCR扩增，配液参见表1/2。

| NEB Q5 Hot Start DNA Polymerase | Pool1 | Pool2 |
| :-------- | :--------:| :--------:|
| 5X Q5 Reaction Buffer | 5 µl | 5 µl |
| 10 mM dNTPs | 0.5 µl | 0.5 µl |
| Q5 Hot Start DNA Polymerase | 0.25 µl | 0.25 µl |
| Primer Pool 1 or 2 (10µM) | 3.6 µl | 3.6 µl |
| Nuclease-free water | 13.15 µl | 13.15 µl |
| cDNA | 2.5uL | 2.5uL |
| **Total** | **25 µl** | **25 µl** |

| Qiagen Multiplex PCR Kit | Pool1 | Pool2 |
| :-------- | :--------:| :--------:|
| 2X Qiagen Multiplex PCR Master Mix | 12.5 µl | 12.5 µl |
| Primer Pool 1 or 2 (10µM) | 3.6 µl | 3.6 µl |
| 5X Q Solution | 2.5 µl | 2.5 µl |
| Nuclease-free water |  3.9 µl | 3.9 µl |
| cDNA | 2.5uL | 2.5uL |
| **Total** | **25 µl** | **25 µl** |

>根据体系优化，每个引物的终浓度为0.015uM。

- NEB Q5 的PCR程序为：

| 序号 | 温度 | 时间 | 前往 |
| :-------- | :--------:| :--------:| :--------:|
| 1. | 98°C | 30s | |
| 2. | 98°C | 15s | |
| 3. | 65°C | 5min | 前往2, 循环25～35cycles |
| 4. | 4°C | Hold | |

- Qiagen Multiplex Kit的PCR程序为：

| 序号 | 温度 | 时间 | 前往 |
| :-------- | :--------:| :--------:| :--------:|
| 1. | 95°C | 15min | |
| 2. | 94°C | 30s | |
| 3. | 58°C | 90s | |
| 4. | 72°C | 45s | 前往2, 循环25～35cycles |
| 5. | 72°C | 10min | |
| 6. | 4°C | Hold | |

>Ct18-21的样本用25循环, Ct 35的样本用35个循环。CT值介于21～35的根据分布自行设置。

>应用本方法到其他病原检测时，需要设计多重PCR的话，可以使用[这里](http://primal.zibraproject.org)的引物设计工具获得引物组。


## 3. Amplicon纯化回收

### 3.1 试剂与耗材

- [Elution Buffer](https://www.qiagen.com/us/products/discovery-and-translational-research/lab-essentials/buffers-reagents/buffer-eb)
- [AMPure XP Beads](https://www.beckman.com/reagents/genomic/cleanup-and-size-selection/pcr)

### 3.2 实验流程

- 将每个样本的Pool1和Pool2混合，因为测序文库与上样量有关，本实验流程设计的是针对8个以上样本的上样量。建议检测样本带质控DNA做为阴性对照，以便检查是否有污染序列。
- 将AMPure XP Beads 提前30min取出置于常温，使用前在震荡仪上混合均匀。
- 取50uL AMPure XP Beads 到 50uL 混合pooling Amplicon中。震荡混匀（增加震荡时间能极大提高回收率），离心甩干。

> 1x 磁珠吸附250～1000bp效率较好。具体[参见](https://research.fhcrc.org/content/dam/stripe/hahn/methods/mol_biol/SPRIselect%20User%20Guide.pdf)

- 室温放置5min。
- 将EP管放置于磁力架上2min，直至液体澄清。
- 吸弃液体，用200uL新鲜配置的70～80%乙醇洗2次。
- 用10uL移液器吸弃残液，开盖1min，让乙醇挥发。
- 从磁力架上取下，用15～30uL Elution Buffer或者水重悬磁珠，用枪头或手指轻轻混匀，静置2min。
- 置于磁力架上，液体澄清后吸取30uL液体到一个新的EP管中。取1uL进行Qubit定量。

>样本浓度与PCR产物回收相关性：
> - 10E5大约200ng/uL
> - 10E4大约120ng/uL
> - 10E3大约80ng/uL
> - 10E2大约25ng/uL
> - 10E1大约10ng/uL

## 4. 文库构建

### 4.1 Nanopore 测序文库

#### 4.1.1 试剂与耗材

- Nanopore EXP-NBD104/114
- Nanopore SQK-LSK109
- Nanopore 测序芯片制备试剂盒 EXP-FLP002
- Eppendorf DNA LoBind Tubes
- [NEBNext® Ultra™ II DNA Library Prep Kit for Illumina E7645](https://international.neb.com/products/e7645-nebnext-ultra-ii-dna-library-prep-kit-for-illumina)

#### 4.1.2 实验流程

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

### 4.2 illumina Miseq 测序文库

#### 4.2.1 试剂与耗材

- Nextera XT LIbrary Preparation Kit
- NEB Ultra II DNA Library Prep Kit for Illumina
- AMPure XP Beads
- Miseq Reagent v2 PE150

> 对于浓度较高的样本，illumina扩增子测序推荐使用PCR-Free的方式进行。这里用NEB Ultra II DNA Library Prep Kit for Illumina(E7645)文库构建试剂盒进行。如果使用Nextera XT之类的酶片段化试剂盒，后期数据分析时要考虑引物扩增位点对于consensus序列的影响。

#### 4.2.2 实验流程

- 将 pool1 和 pool2 的 PCR 产物混合后，使用Qubit进行定量。NEB文库构建试剂盒适合500pg~1ug的起始input DNA量。
- 取新的 PCR 管，按照下表配液。用200uL枪的50uL量程吹吸10次混匀。离心甩干。

| 试剂 | 量 |
| :----- | :--: |
| NEBNext Ultra II End Prep Enzyme Mix | 3 µl |
| NEBNext Ultra II End Prep Reaction Buffer | 7 µl |
| DNA | 50 µl |
| **总量** | **60 µl** |

- 放置于PCR仪器上，热盖温度大于等于75°C ，运行以下程序：20°C 30min，65°C 30min，4°C hold。
- 按照下表稀释 adaptor。建议初始DNA量大于100ng，如200ng左右。

| Input DNA | Adaptor | Adaptor 工作浓度 |
| :----- | :--: | :--: |
| 1ug-101ng | 不稀释 | 15uM |
| 100ng-5ng | 1:10稀释 | 1.5 uM |
| <5ng | 1:25稀释 | 0.6uM |

- 根据下表配置Adaptor连接体系。NEBNext Ultra II Ligation Master Mix加入前要颠倒混匀。

| 试剂 | 量 |
| :----- | :--: |
| 前一步反应液 | 60 µl |
| NEBNext Ultra II Ligation Master Mix | 30 µl |
| NEBNext LIgation Enhancer | 1 µl |
| NEBNext Adaptor for illumina | 2.5 µl |
| **总量** | **93.5 µl** |

>可以采用的Adaptor产品：单端（E7350），双端（E7335,E7500,E7710,E7730,E7600,E7535,E6609）

- 用200uL枪的80uL量程吹戏10次混匀。离心甩干。
- 置于金属浴或PCR仪上（不开热盖）20°C 15min。
- 向反应体系中加入3uL USER Enzyme（这个试剂包含在接头试剂盒中）。在PCR仪器上37°C 15min（热盖大于等于47°C）。
- 将液体转移到深孔板内，当input DNA大于50ng时，进行下面的DNA纯化回收。对于小于50ng DNA的样本，只需要用磁株纯化1次，用87uL量。
- 加入18uL AMPure Beads到反应液中，吹吸10次混匀，注意枪头中的液体要全部打出。室温孵育5min。
- 将深孔板置于磁力架上，静置5min待液体澄清。吸出液体到新的孔中。
- 将深孔板从磁力架上取下，加入10uL AMPure Beads，吹吸10次混匀后，室温孵育5min。
- 将深孔板置于磁力架上，静置5min待液体澄清。将液体吸弃。
- 加入80%新鲜配置的乙醇200uL，室温放置30s后吸弃。重复洗1次。
- 开盖风干磁珠5min。如果风干时间过长，会导致DNA回收率下降。
- 从磁力架上取下深孔板，17uL 10mM Tris-HCl或0.1xTE洗脱。
- 震荡仪上1800rpm震荡10min，室温放置2min。
- 将深孔板置于磁力架上，5min后待液体澄清，取1uL进行Qubit定量。
- 吸取15uL液体到新的PCR板。


## 5. 数据分析

### 5.1 Nanopore 测序数据

#### 5.1.1 软件安装

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

#### 5.1.2 分析流程

> 与illumina原始数据是图像文件不同，Nanopore的原始数据以HDF5格式保存。HDF是"Hierarchical Data Format"首字母缩写，HDF5是一种纯文本，类似json的数据格式记录电信号。查看HDF5原始数据可以用 [hdf5_tools](https://support.hdfgroup.org/products/hdf5_tools/) 工具。数据文件后缀一般用.fast5表示。

Nanopore 测序数据期初由于其错误率高而“闻名”，特别是对二聚体的区分。电信号数据随着处理软件和算法的不断改进，准确率得到不断的提升。HDF5数据进行basecalling的软件目前很多，采用的算法也各异。官方的如albacore，guppy等，第三方工具也有很多。

MinKnow软件 将fast5数据复制到服务器，使用guppy进行分析

guppy是一个是用

我们的实验室用一台Workstation挂载MinION测序仪，进行测序工作。基本配置是i7-7700k+16G+1TSSD/1T HD。一般可以做到70%以上的MinKnow实时basecalling。但是如果要使用guppy进行更准确的basecalling的话，就需要将数据同步到服务器上进行。

```bash
# 在服务器上 rsync 同步数据到本地
(ncov)$ rsync -avz username@minknow_ip:/var/lib/minknown/path/to/fast5 .
# 或者在工作站上同步数据到服务器端
$ rsync -avz /var/lib/minknown/path/to/fast5 username@server_ip:~/data/fast5
```

这里设置服务器端 fast5 数据目录存放于：$HOME/data/fast5；采用的是R9.4的测序芯片，各个不同测序芯片basecalling设置参数可以参见[这里](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html)

```bash
(ncov)$ guppy_basecaller -c dna_r9.4.1_450bps_fast.cfg \
> -i $HOME/data/fast5 -s run_name -x auto -r
```

>我们的MinION测序仪连接的是一台Ubuntu 16.04的Workstation，测序仪最大测序状态时，跑默认的basecalling基本上可以实时达到80%的状态。如果需要用guppy做更准确的basecalling，我们会用inotify-tools之类的工具加上 rsync 实时同步到服务器，在服务器上进行。


artic 的 demultiplex 命令是

```bash
$ porechop --verbosity 2 --untrimmed -i "re_all.fastq" -b ./tmp5vx2iwn0 --native_barcodes --discard_middle --require_two_barcodes --barcode_threshold 80 --threads 8 --check_reads 10000 --barcode_diff 5 > re_all.fastq.demultiplexreport.txt
```


### 5.2 Illumina 测序数据

#### 5.2.1 软件安装

>暴发疫情中病毒毒株发生结构变异的可能性较小，可以直接使用比对参考基因组的方法获得毒株基因组序列。如果是研究新发的远源病毒或长时间演化的病毒序列，可能与参考基因组差异较大时，应结合组装与比对的方法获得基因组序列。

```bash
$ conda create -n mapping
$ conda activate mapping
(mapping)$ conda install bwa samtools igv
```

#### 5.2.2 分析流程

- 使用Miseq自带的SAV软件进行basecalling，生成fastq数据复制到服务器端
- 可以采用bwa+samtools进行参考序列比对获得一致性序列，或者采用bwa等工具比对获得新冠病毒序列后在进行de novo组装。

```bash
# 最基本的比对分析流程
(mapping)$ bwa index MN909847.3.fasta
(mapping)$ bwa mem -k 45 -R "@RG\tID:HZ-1_1\tSM:HZ-1" MN909847.3.fasta \
> HZ-1_S1...R1.fastq.gz HZ-1_S1...R2.fastq.gz | samtools view -bS > HZ-1.bam
(mapping)$ samtools sort -O bam -o HZ-1.sorted.bam -@4 HZ-1.bam
(mapping)$ samtools index HZ-1.sorted.bam
(mapping)$ igv HZ-1.sorted.bam

# 以bwa为例进行reads过滤后spades组装
$ bwa index MN908947.3.fasta
$ bwa mem -k 47 -R "@RG\tID:HZ-1\tSM:HZ-1" MN908947.3.fasta \
> HZ-1_S1_L001_R1_001.fastq.gz HZ-1_S1_L001_R2_001.fastq.gz -t 40 | \
> samtools view -bS - | samtools sorted -@40 -O bam HZ-1.sorted.bam
$ samtools index HZ-1.sorted.bam
$ samtools faidx MN908947.3.fasta
$ samtools mpileup -uf MN908047.3.fasta HZ-1.sorted.bam | \
> bcftools call --ploidy 1 -mv -Oz > HZ-1.vcf.gz
$ tabix HZ-1.vcf.gz
$ cat MN908947.3.fasta | bcftools consensus HZ-1.vcf.gz > Hz-1_consensus.fasta

# de novo 组装
$ bwa index MN908947.3.fasta
$ bwa mem -k 45 -R "@RG\tID:HZ-1\tSM:HZ-1" MN908947.3.fasta \
> HZ-1_S1_L001_R1_001.fastq.gz HZ-1_S1_L001_R2_001.fastq.gz -t 4 | \
> samtools view -bS > HZ-1.bam
# 提取双端都比对到参考序列的reads
$ samtools view -bF 12 HZ-1.bam > HZ-1.mapped.bam
# bam格式转化成fastq格式
$ bamToFastq -i HZ-1.mapped.bam -fq HZ-1.mapped_R1.fastq -fq2 HZ-1.mapped_R2.fastq
# de novo assembly
$ shovill --trim --outdir test --R1 HZ-1.mapped_R1.fastq -fq2 HZ-1.mapped_R2.fastq --ram 16 --cpus 40
$ bwa mem -t 40 -x ont2d MN908947.3.fasta test/contigs.fa | samtools view -bS - | samtools sort -o test.sorted.bam -
$ samtools index test.sorted.bam
```

## 自定义分析流程

针对MinKnow生成的数据

```bash
# 生成样本的 fastq
$ for i in *.fastq; do zcat $i >> s1.fq; done
$ gzip s1.fq
# 去除可能的接头序列
$ porechop -t 4 -i s1.fq.gz --format fastq.gz -o s1.clean.fq.gz
# 去除低质量序列
$ gunzip -c s1.clean.fq.gz | NanoFilt -q 10 -l 400 --maxlength 700 | gzip > s1.clean.hq.fq.gz
# 获得比对结果
$ bwa index MN908947.3.fasta
$ bwa mem -t 4 -x ont2d MN908947.3.fasta s1.clean.hq.fq.gz | \
> samtools view -bS - | \
> samtools sort -o s1.clean.hq.sorted.bam -
$ samtools index s1.clean.hq.sorted.bam
$ samtools faidx MN908947.3.fasta
# vcf calling
$ bcftools mpileup -f --threads 4 MN908047.3.fasta s1.clean.hq.sorted.bam | \
> bcftools call --ploidy 1 -mv -Oz > s1.vcf.gz
$ tabix s1.vcf.gz
# 获得 consensus 序列
$ bcftools consensus -f MN908947.3.fasta -H 1 s1.vcf.gz -o s1_consensus.fasta


# 或者使用minimap工具进行比对

```

##  参考资料

- https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w
- https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
- https://nanoporetech.com/resource-centre/nanopore-sequencing-sars-cov-2-genome-introduction-protocol
