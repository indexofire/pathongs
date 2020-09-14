# bwa + samtools 基本分析流程

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

可以采用基于比对的方法获得参考基因组 consensus 序列。这种方法可以应用在病毒等基因组序列的获得。

- 比对软件: bwa
- BAM软件: samtools
- 去重复: gatk4 MarkDuplicates
- 生成vcf: freebayes
- 注释vcf: snpEff

## 命令操作

```bash
# 对 fasta 格式的参考序列 ref.fa 建立索引
$ bwa index ref.fa

# 利用 mem 算法将 illumina reads 比对到参考基因组上
# 并使用 samtools 将其转换为二进制 bam 个是文件
$ bwa mem -t 40 ref.fa -R "@RG\tID:1\tSM:S1" 1_S1_L001_R1.fastq.gz 1_S1_L001_R2.fastq.gz | \
> samtools view -bS > 1.bam

# 将 bam 文件排序
$ samtools sort -O bam 1.sorted.bam 1.bam

# 使用 gatk 去除重复需里
$ gatk MarkDuplicates -I 1.sorted.bam -O 1.markdup.sorted.bam -M 1.markdup.sorted_metrics.txt

# 对排序并去重复的序列建立索引
$ samtools index 1.markdup.sorted.bam

# 使用 freebayes 将 bam 格式的结果转换为 vcf 格式的文件
$ freebayes -m 20 -p 1 -f ref.fa 1.markdup.sorted.bam > 1.vcf

# 使用 snpEff 获得 SNP 为点注释
$ snpEff lmo_database 1.vcf > 1.snpEff.vcf
```


## Reference

[BWA中文手册](http://cncbi.github.io/BWA-Manual-CN/)
