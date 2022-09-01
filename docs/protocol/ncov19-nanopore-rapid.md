# 新冠病毒快速法纳米孔测序

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! abstract "内容简介"
    本实验流程是根据纳米孔测序方法中快速法来对 Sars-Cov-2 基因组纳米孔测序的方案。同时根据我们实验室实际工作中应对疫情和监测工作的经验生成的实验方案。本实验流程适合针对12个样本以下的检测工作并希望快速获得结果的应用，比较适合应急疫情的处置。

## 1. 实验方案

对于单个样本，采用SQK-RAD004试剂盒，对于2～12个样本采用SQK-RBK004试剂盒构建文库。gDNA获得的方式采用AN、VS和MD进行PCR扩增，然后混合扩增子并纯化回收后得到。对于常见的奥密克戎变异株，我们常采用AN+MD 2组引物进行扩增，增加覆盖度。

!!! note "新版试剂"
    Oxford Nanopore公司新出的RBK112试剂盒待我们到货试用后再进行确认。

## 2. 实验流程

### 2.1 cDNA获取

参见low-cost，由于产物进行纯化回收，因此PCR扩增的循环数可以根据CT值适当调整。

### 2.2 文库构建

#### 2.2.1 RAD004

取7.5uL（浓度约为50～60ng/uL）纯化定量后的扩增子，加入2.5uL FRA，混合后置于PCR仪上30°C 1min，80°C 1min。然后立即置于冰上至少1min。加入1uL RAP 混合后20°C放置5min。按下表配液并上机。

| 试剂 | 体积 |
| :-- | :--:|
| SQB | 34 µL |
| LB | 25.5 µL |
| Nuclease-free water | 4 µl |
| 文库 | 11 uL |
| **Total** | **75 µl** |

#### 2.2.2 RBK004

2～12个样本每个扩增子取7.5uL(浓度约为15～30ng/uL) 纯化定量后的扩增子，加入2.5uL RBXX(XX为01-12，12个样本对应12个barcodes)，混合后置于PCR仪上30°C 1min，80°C 1min。然后立即置于冰上至少1min。

将所有文库pooling在一个lo-binding管中，加入0.8x磁珠（8uL~96uL），混匀后室温放置5min，置于磁力架上2min，吸弃上清液，用80%乙醇洗2次。用10uL移液器吸弃残液，开盖风干1min。加入10uL Nuclease-free 水，室温放置2min。置于磁力架上2min后吸取上清液至PCR管中，加入1uL RAP 混合后20°C放置5min。按RARD004配液表上机。

## 3. 数据分析

分析流程我们采用自己编写的流程脚本进行，参见[这里]

```bash
# 对于RAD004文库: run_rad004 samplename sampledate
(artic)$ run_rad004 22cov100 2022-02-10
# 结果序列文件在上一及目录 22cov100.fasta

# 对于RBK004文库，需要建立barcodes.csv
# barcodes.csv
barcode,sample,lab,date
barcode01,22cov100,HZCDC,2022-01-05
barcode02,22cov101,HZCDC,2022-01-06

# 运行分析流程
(artic)$ run_rbk004
# 结果序列文件在上一及目录
(artic)$ ls ../
22cov100.fasta 22cov101.fasta
```

**run_rad004**

```bash
#!/usr/bin/env bash
# Usage
# $ run_rad004 22cov1500 2022-03-05
# get genomic fasta in the upper level directory, 22cov1500.fasta
#

ls fastq_pass/*.gz | parallel -k -j8 zcat {} > $1.fastq
cat $1.fastq | NanoFilt -q 10 --maxlength 1200 --headcrop 25 --tailcrop 25 > $1.clean.fastq
minimap2 -ax map-ont -t 8 -R '@RG\tID:'$1'\tSM:'$1 ref.fasta $1.clean.fastq | samtools view -bS -F 4 - | samtools sort -o $1.sorted.bam
samtools index $1.sorted.bam
medaka consensus --model r941_min_sup_g507 --chunk_len 1400 --chunk_ovlp 200 $1.sorted.bam $1.hdf
medaka variant ref.fasta $1.hdf $1.vcf
bgzip -f $1.vcf
tabix -f -p vcf $1.vcf.gz
longshot -P 0 -F -A --no_haps --bam $1.sorted.bam --ref ref.fasta --out $1.longshot.vcf -v $1.vcf.gz
artic_vcf_filter --medaka $1.longshot.vcf $1.pass.vcf $1.fail.vcf
bgzip -f $1.pass.vcf
tabix -p vcf $1.pass.vcf.gz
artic_make_depth_mask --store-rg-depths ref.fasta $1.sorted.bam $1.coverage_mask.txt
artic_mask ref.fasta $1.coverage_mask.txt $1.fail.vcf $1.preconsensus.fasta
bcftools consensus -f $1.preconsensus.fasta $1.pass.vcf.gz -m $1.coverage_mask.txt -o $1.consensus.fasta
echo '>SARS-CoV-2/human/CHN/HZCDC-'$1'/'$2 > ../$1.fasta
sed -n '2,$p' $1.consensus.fasta >> ../$1.fasta
```

**run_rbk004**

```bash
#!/usr/bin/env bash
# Usage
# $ run_rbk004
# get genomic fasta in the upper level directory which we used to put in
#

REF=$CONDA_PREFIX/ref.fasta

if [ ! -f barcodes.csv ];then
    echo "No barcodes.csv file found"
else
    while IFS=, read -r barcode sample lab date
    do
    if [ $barcode = "barcode" ]
    then
        continue
    else
    ls fastq_pass/${barcode}/*.gz | parallel -k -j8 zcat {} > ${barcode}.fastq
    porechop -t 8 -i ${barcode}.fastq | NanoFilt -q 10 --maxlength 1200 --headcrop 25 --tailcrop 25 > ${barcode}.clean.fastq
    minimap2 -ax map-ont -t 8 -R '@RG\tID:'$sample'\tSM:'$sample $REF ${barcode}.clean.fastq | samtools view -bS -F 4 - | samtools sort -o ${barcode}.sorted.bam
    samtools index ${barcode}.sorted.bam
    medaka consensus --model r941_min_sup_g507 --chunk_len 1400 --chunk_ovlp 200 ${barcode}.sorted.bam ${barcode}.hdf
    medaka variant $REF ${barcode}.hdf ${barcode}.vcf
    bgzip -f ${barcode}.vcf
    tabix -f -p vcf ${barcode}.vcf.gz
    longshot -P 0 -F -A --no_haps -b ${barcode}.sorted.bam -f $REF -o ${barcode}.longshot.vcf -v ${barcode}.vcf.gz
    artic_vcf_filter --medaka ${barcode}.longshot.vcf ${barcode}.pass.vcf ${barcode}.fail.vcf
    bgzip -f ${barcode}.pass.vcf
    tabix -p vcf ${barcode}.pass.vcf.gz
    artic_make_depth_mask --store-rg-depths ref.fasta ${barcode}.sorted.bam ${barcode}.coverage_mask.txt
    artic_mask $REF ${barcode}.coverage_mask.txt ${barcode}.fail.vcf ${barcode}.preconsensus.fasta
    bcftools consensus -f ${barcode}.preconsensus.fasta ${barcode}.pass.vcf.gz -m ${barcode}.coverage_mask.txt -o ${barcode}.consensus.fasta
    echo '>SARS-CoV-2/human/CHN/'${lab}'-'${sample}'/'${date} > ../${sample}.fasta
    sed -n '2,$p' ${barcode}.consensus.fasta >> ../${sample}.fasta
    fi  
    done < barcodes.csv
fi
```
