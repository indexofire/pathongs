# GATK4

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 安装

### 1.采用 conda 安装

```bash
# 创建虚拟环境
$ conda create -n gatk4
$ conda activate gatk4

# 安装 gatk4
(gatk4)$ conda install gatk4

# 显示 gatk4 工具列表
(gatk4)$ gatk --list

# 查看二级命令的帮助内容，这里 TOOLS 为各个二级命令名称
(gatk4)$ gatk TOOLS --help
```

### 2.采用 docker 安装

这里使用官方镜像地址: https://hub.docker.com/r/broadinstitute/gatk/ 。

1. 首先在 [docker hub](https://hub.docker.com/signup) 注册一个账户，
2. 在终端用`docker login`命令登录

```bash
# 使用注册的用户名
$ docker login -u your_username -p your_password
$ docker pull broadinstitute/gatk

# 以gatk为镜像，启动容器并进入容器的交互界面
$ docker run -it -rm broadinstitute/gatk bash
root@bb3d14c7da19:/gatk#

# 看到root@bb3d14c7da19:/gatk# 表示成功进入容器
# 参数含义
-i 交互界面操作
-t 终端中执行
--rm 表示推出容器后立即删除该容器，避免占据磁盘空间

# 查看二级命令的帮助内容，这里 TOOLS 为各个二级命令名称
root@bb3d14c7da19:/gatk# gatk TOOLS --help
```

```bash
# 如果你的docker没有启动，在Archlinux下可以systemctl启动
$ sudo systemctl start docker.service

# 如果出现该错误，因为使用了新的内核，重新启动即可
# failed to create endpoint on network bridge: failed to add the host sandbox pair interfaces operation not supported.
```

## 使用

### 1.GATK 工具

查看工具列表

```bash
$ gatk --list
```

GATK4的工具列表分类为：

- Base Calling
    * CheckIlluminaDirectory:
    * CollectIlluminaBasecallingMetrics
    * CollectIlluminaLaneMetrics
    * ExtractIlluminaBarcodes
    * IlluminaBasecallsToFastq
    * IlluminaBasecallsToSam
    * MarkIlluminaAdapters
- Copy Number Variant Discovery
- Coverage Analysis
- Diagnostics and Quality Control
- Genotyping Arrays Manipulation
- Intervals Manipulation
- Metagenomics
- Methylation-Specific Tools
- Other
- Read Data Manipulation
- Reference
- Short Variant Discovery
- Structural Variant Discovery
- Variant Evaluation and Refinement
- Variant Filtering
    * FixVcfHeader
    * GatherVcfs
    * GatherVcfsCloud
    * LeftAlignAndTrimVariants
    * LiftoverVcf
    * MakeSitesOnlyVcf
    * MakeVcfSampleNameMap
    * MergeVcfs
    * PrintVariantsSpark
    * RemoveNearbyIndels
    * RenameSampleInVcf
    * SelectVariants
    * SortVcf
    * SplitVcfs
    * UpdateVCFSequenceDictionary
    * UpdateVcfSequenceDictionary
    * VariantAnnotator
    * VcfFormatConverter
    * VcfToIntervalList

### 2.突变位点分析流程

这里以副溶血性弧菌为例，参考基因组为RMID2210633

```bash
# 弧菌有2个环状染色体，这里将 RMID 2210633 的2条染色体序列分别构建索引
# 也可以单独一个fasta格式基因组序列构建索引
$ samtools faidx NC_004603.fna
$ samtools faidx NC_004605.fna

# 用bwa参考基因组构建索引
$ bwa index NC_004603.fna
$ bwa mem -t 4 -R '@RG\tID:test\tSM:RMID2210633' NC_004603.fna \
> SRR8535443_1.fastq.gz SRR8535443_2.fastq.gz | samtools view -bS > SRR8535443_chr1.bam
$ samtools sort -@4 -o SRR8535443_chr1_sorted.bam SRR8535443_chr1.bam

# 去除光学重复序列
$ gatk MarkDuplicates -I SRR8535443_chr1_sorted.bam \
> -O SRR8535443_chr1_sorted_markdup.bam \
> -M SRR8535443_chr1_sorted_markdup-metrics.txt

# bam文件建立索引
$ samtool index SRR8535443_chr1_sorted_markdup.bam

# gatk 对参考序列生成 .dict 文件
$ gatk CreateSequenceDictionary -R NC_004603.fna

$ gatk HaplotypeCaller -R NC_004603.fna --emit-ref-confidence GVCF \
> -I SRR8535443_chr1_sorted_markdup.bam \
> -O SRR8535443_chr1.gvcf

$ gatk GenotypeGVCFs -R NC_004603.fna \
> -V SRR8535443_chr1.gvcf
> -O SRR8535443_chr1.vcf
```

!!! Note "注意事项"
    gatk4 默认流程call的 SNP 数量比`freebayes`、`snippy`要更多。主要原因是freebayes把许多4个长度碱基1/4位置变化的序列作为重组complex而不是snp，而gatk则作为snp处理。另外gatk使用的reads数跟多，可以call出头尾和一些低覆盖度的区域。但是gatk会有一些无覆盖区域的snps出现，需要进一步研究是什么原因造成的。如果无法消除，则gatk可能不适合作为微生物的snp calling工具。
