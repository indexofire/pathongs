# GATK3

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 突变位点分析流程

```bash
$ bwa index ref.fa
$ bwa mem -t 4 -R '@RG\tS1SM\tS1' ref.fa S1_R1.fq.gz S1_R2.fq.gz > S1.sam
$ samtools view -S -b S1.sam > S1.bam
$ samtools sort -@ 4 -O bam -o S1.sorted.bam S1.bam
$ java -jar picard.jar MarkDuplicates \
> I=S1.sorted.bam \
> O=S1.sorted.markdup.bam \
> M=S1.markdup_metrics.txt
$ samtools index S1.sorted.markdup.bam
```
