# 分析流程1

- 比对软件: bwa
- BAM软件: samtools
- 去重复: gatk4 MarkDuplicates
- 生成vcf: freebayes
- 注释vcf: snpEff 

```bash
$ bwa index ref.fa
$ bwa mem -t 40 ref.fa -R "@RG\tID:1\tSM:S1" 1_S1_L001_R1.fastq.gz 1_S1_L001_R2.fastq.gz | \
> samtools view -bS > 1.bam
$ samtools sort -O bam 1.sorted.bam 1.bam
$ gatk MarkDuplicates -I 1.sorted.bam -O 1.markdup.sorted.bam -M 1.markdup.sorted_metrics.txt
$ samtools index 1.markdup.sorted.bam


$ freebayes -m 20 -p 1 -f ref.fa 1.markdup.sorted.bam > 1.vcf
$ snpEff lmo_database 1.vcf > 1.snpEff.vcf


```



