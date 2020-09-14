# 病毒基因组组装

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

与细菌一般可以培养从而获得高质量的基因组DNA，然后进行基因组测序不同。有些病毒难以培养（如诺如病毒），有些病毒培养生物安全等级要求较高，且时间长。
病毒基因组测序有很多方法，这里只针对二代高通量测序数据进行讨论。用二代测病毒基因组有不同的方法，比如对临床样本中已知病毒进行特异性PCR扩增后测序；有些则采用随机引物扩增后再测序，这种方法可以针对未知病毒序列，但由于不同样本（比如临床样本）中含有的干扰核酸较多（如人源DNA等），往往需要实验中进行特殊处理或者生物信息学分析时筛选。

因此我们只讨论特异性PCR扩增后（类似一代测序）用二代测序分析已知病毒的方法。和细菌单克隆培养物不同，临床样本即使PCR扩增后，

由于病毒基因组位点变异较明显，但是大片段的异位重组比较少，因此可以采用基于 mapping 的方法获得基因组。

首先要选择参考序列，由于病毒编译大，我们一般尽可能选择近源的毒株。我么那可以采用 denovo 的方法先获得一个大片段高覆盖的非污染序列，然后 blast 找到近源毒株基因组序列，将其作为参考基因组再使用 mapping 的方法获得测序毒株的基因组。

## 采用 freebayes

```bash
$ spades.py
$ blast
$ bwa index ref.fasta
$ bwa mem
$ samtools index
$ samtools
$ freebayes
```

## 采用 GATK4

```bash
$ trimmomatic

# ID:
# PL: Platform 测序平台，这里数据来自 illumina 平台
# PU:
# LB:
# SM: Sample 样本名称
$ bwa mem -t 40 -M -Y -R "@RG\tID:S1\tPL:illumina\tPU:\tLB:\tSM:" ref.fasta \
>
$ samtools sort -@ 4 -m 16G -O bam -o .sorted.bam

$ gatk MarkDuplicates -I -M -O
$ samtools index sorted.markdup.bam
$ gatk BaseRecalibrator -R ref.fasta -I sorted.markdup.bam \
> -O sorted.markdup.recal.table
$ gatk ApplyBQSR --bqsr-recal-file sorted.markdup.recal.table -R ref.fasta
> -I sorted.markdup.bam -O sorted.markdup.BQSR.bam
$ samtools index
$ gatk HaplotypeCaller -R ref.fasta -I sorted.markdup.BQSR.bam \
> -O HC.vcf.gz

# SNP call
$ gatk VariantRecalibrator
# Indel call

```

## 使用 snippy
