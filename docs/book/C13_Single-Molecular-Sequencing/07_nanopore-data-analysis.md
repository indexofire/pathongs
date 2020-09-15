# Nanopore 数据处理

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

以针对一个1.5k左右的病毒扩增子测序文库进行的 Nanopore 测序数据处理。


```
# 目录结构
fastq_pass
├── barcode01
│   ├── 00_raw              # nanopore guppy basecalling 获得的fastq_pass数据
│   ├── 01_centrifuge       # 用centrifuge对数据进行扫描后获得的基本物种分布
│   ├── 02_porechop         # 用porechop去除barcodes的reads
│   ├── 03_nanofilt         # 用nanofilt去除低质量序列和短序列        
├── barcode02           
├── ...
└── barcode24
```

## 测序结果扫描

```bash
(centrifuge)$ export DB=$HOME/dbs/centrifuge/b_compressed+h+v
(centrifuge)$ ls 00_raw/*.fastq | parallel --max-args=1 centrifuge -x $DB -U {1} -S 01_centrifuge/{1}.result --report-file 01+centrifuge/{1}.reprot
```

## 去除接头

使用 porechop 去除 barcodes，如果在 MinKnow 中使用 guppy 去除过一次，建议再操作一次。

```bash
(nanopore)$ ls 00_raw/*.fastq | parallel --max-args=1 porechop -i {1} -o 02_porechop/{1}
```

## 去除低质量序列和短片段

对于 Nanopore 生成的碱基质量做一便过滤，这里使用Q值为9作为合格线，reads 长度超过1200的保留。

```bash
(nanopore)$ ls 02_porechop/*.fastq | parallel --max-args=1 NanoFilt -q 9 -l 1200 {1} | gzip > 03_nanofilt/trimed.fastq.gz
```

## 去除引物区域

## 直接组装获得基因组


## 比对到参考基因组获得基因组

采用 minimap2

```bash
(nanopore)$ minimap2
```

采用 bwa mem

```bash
(nanopore)$ bwa index ref.fa
(nanopore)$ bwa mem -t 4 -R "@RG:S1\tSM:S1" ref.fa 03_nanofilt/trimed.fastq.gz | \
> samtools view -bS - | samtools sort -@4 -O bam -o sorted.fastq.gz
(nanopore)$ bwa index sorted.fastq.gz
```


## 获得基因组及突变位点

## 将整个程序构建 snakemake 流程

## Reference
