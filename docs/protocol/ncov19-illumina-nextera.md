# illumina Miseq 平台进行 SARS-CoV-2 病毒基因组测序

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! abstract "内容简介"
    简介实验室采用二代测序方法在 illumia miseq 平台进行 Sars-Cov-2 基因组测序方案，结合ivar分析流程获得基因组一致性序列。

## 实验方法

- 片段化法： Nextera XT 文库试剂盒
- 连接法： NEBNext UltraII 文库试剂盒

2种方法前期逆转录和PCR扩增基本一致。对于片段化法建议使用MD或者VS long引物进行扩增实验。

### 1. 逆转录

SuperScript IV 试剂盒进行逆转录。逆转录产物cDNA可在4°C过夜或者-20°C放置一周。

### 2. PCR扩增

- Pool1 + Pool2 混合后，加入0.8x AMPure Beads，吹吸混匀，室温放置10min。
- 置于磁力架上，5min。移液器吸弃残液。
- 加入180uL 80%乙醇，30s后吸弃。共洗2次，不要将离心管从磁力架上取下。
- 将离心管取下，离心甩干，放回磁力架上，将残液吸弃，风干1min。
- 取下离心管，加入30uL水，重悬磁珠，室温放置5min。
- 置于磁力架上2min，将液体吸出即为纯化回收的扩增DNA。
- 取1uL进行定量。

### 3. 文库构建

#### 3.1 连接法

**试剂与耗材**

- NEB Ultra II DNA Library Prep Kit for Illumina

**实验步骤**

1. 根据定量的扩增子，计算所要进行文库构建的input DNA量，保证每个样本大于100ng，体积为26uL，相当与至少浓度为4ng/uL。

| Reagent | Volume |
| ------- | ------ |
| NEBNext Ultra II FS Reaction Buffer | 7uL |
| NEBNext Ultra II FS Enzyme Mix | 2uL |


#### 3.2 片段化法

**试剂与耗材**

- Nextera XT LIbrary Preparation Kit
- AMPure XP Beads
- Miseq Reagent v2 PE150

**实验步骤**

1. 将 pool1 和 pool2 的 PCR 产物混合后，使用Qubit进行定量。NEB文库构建试剂盒适合500pg~1ug的起始input DNA量。
2. 取新的 PCR 管，按照下表配液，混匀后离心甩干。

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

**数据分析**

>暴发疫情中病毒毒株发生结构变异的可能性较小，可以直接使用比对参考基因组的方法获得毒株基因组序列。如果是研究新发的远源病毒或长时间演化的病毒序列，可能与参考基因组差异较大时，应结合组装与比对的方法获得基因组序列。

- 使用Miseq自带的SAV软件进行basecalling，生成fastq数据后复制到安装分析软件的服务器端
- 可以采用bwa+samtools进行参考序列比对获得一致性序列，或者采用bwa等工具比对获得新冠病毒序列后在进行de novo组装。

```bash
# 最基本的比对分析流程
$ bwa index MN909847.3.fasta
# 去除可能含有的引物序列
$ fastp Sample1_R1.fastq.gz
# 比对到参考基因组
$ bwa mem -R "@RG\tID:seq1\tSM:sample1" MN909847.3.fasta \
> Sample1_R1.fastq.gz Sample1_R2.fastq.gz | samtools view -bS | \
> samtools sort -O bam -o Sample1.sorted.bam -@4
$ samtools index Sample1.sorted.bam
# 打开 igv 可视化查看突变位点
$ igv Sample1.sorted.bam

# 以bwa为例进行reads过滤后spades组装
$ bwa index MN908947.3.fasta
$ bwa mem -R "@RG\tID:seq1\tSM:sample1" MN908947.3.fasta \
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

实验方案采用的是对1200bp的扩增子片段化并构建文库，以PE150模式进行测序。24个样本的逆转录大约需要1h，PCR扩增需要4h，文库构建大约需要12h，Miseq上机测序时间大约26h，下机生成fastq数据以及下游分析大约3h。完成整个过程大约2days时间。

这个流程优点是比较简单稳定，数据质量高，分析过程较为简单，有成熟的成套试剂盒。缺点是测序周期长，如果要加快速度需要采用更快速的机型以及PE75这样的更短的读长进行测序。总体而言更适合回溯研究，而不利于应急疫情快速获得结果。
