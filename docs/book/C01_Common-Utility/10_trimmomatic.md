# 二代数据的trim操作

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## trimmomatic

[trim and filt](https://datacarpentry.org/wrangling-genomics/03-trimming/index.html)

软件版本: v

```bash
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
> SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
> SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
> ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
```
