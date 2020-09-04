# 副溶血性弧菌数据分析

---

!!! Abstract "内容简介"
    本节简述对副溶血性弧菌进行分析的一个项目过程。

## 数据下载

!!! note "软件列表"
    entrez-direct
    sra-tools
    parallel-fastq-dump

选择下载 illlumina 平台的数据，主要针对 miseq 和 hiseq（hiseq 各个平台包括 x ten）。novaseq 目前暂时没有微生物数据。

```bash
# 搜索 SRA 数据库中基因组测序数据的记录，保存成 vpa.sra.xml 文件
$ esearch -query 'txid670[organism] AND "strategy wgs"[properties] AND "biomol dna"[properties]' -db sra | efetch -format docsum > vpa.sra.xml

# 获得 vpa.sra.xml 数据中各个测序平台的数据
$ xtract -input vpa.sra.xml -pattern DocumentSummary -element Platform@instrument_model | sort | uniq -c
      8 454 GS FLX+
      1 454 GS Junior
     18 AB SOLiD 4 System
     25 HiSeq X Ten
    158 Illumina HiSeq 2000
    161 Illumina HiSeq 2500
    991 Illumina MiSeq
    119 Ion Torrent PGM
      2 NextSeq 500
     28 PacBio RS
      1 PacBio RS II
      1 Sequel

# 分离不同测序平台的数据表
# 对于 MiSeq 平台，筛选测序碱基数据大于150M，且没有一个样本对应对多个测序数据的实验
$ xtract -input vpa.sra.xml -pattern DocumentSummary -element \
> Platform@instrument_model Statistics@total_bases Biosample Run@acc | \
> grep MiSeq | awk -F'\t' '{if($2>150000000 && !$5) print $0}' > vpa.sra.miseq.acc
# 确认是否有多个 Run 对应相同的 Biosample
$ awk -F'\t' '{print $3}' vpa.sra.miseq.acc | awk 'a[$0]++'
# 下载获得的 Run accession
$ for i in $(cat vpa.sra.miseq.acc | awk -F'\t' '{print $4}'); do prefetch $i; done
# 提取 biosample 信息
$ for i in $(cat vpa.sra.miseq.acc | awk -F'\t' '{print $3}' | awk '!a[$0]++'); do efetch -db Biosample -id $i > biosample/miseq/$i.xml; sleep 5; done

# 对于 HiSeq 2000, 2500, X Ten 平台的数据，碱基数据大于150M，平均读长超过 PE100 的测序数据
$ cat vpa.sra.xml | xtract -pattern DocumentSummary -element \
> Platform@instrument_model Statistics@total_spots Statistics@total_bases \
> Biosample Run@acc | grep HiSeq | awk -F'\t' '{if($3 >150000000 && $3/$2 > 200 && !$6) print $0}' > vpa.sra.hiseq.acc
# 确认是否有多个 Run 对应相同的 Biosample
$ awk -F'\t' '{print $4}' vpa.sra.hiseq.acc  | awk 'a[$0]++'
# 下载获得的 Run accession
$ for i in $(cat vpa.sra.hiseq.acc | awk -F'\t' '{print $5}'); do prefetch $i; done
# 提取 biosample 信息
$ for i in $(cat vpa.sra.hiseq.acc | awk -F'\t' '{print $4}' | awk '!a[$0]++'); do efetch -db Biosample -id $i > biosample/hiseq/$i.xml; sleep 5; done

# fastq-dump 工具只支持单线程转换，速度收到限制。
# 对于批量 sra 数据，可以用 parallel 来批量进行 fastq 转换
$ parallel "fastq-dump --gzip --split-3" ::: *.sra
# 对于单个 sra 数据，可以利用 parallel-fastq-dump 进行多核加速处理。
$ conda install parallel-fastq-dump
$ parallel-fastq-dump --sra-id ... --threads 40 --split-files --gzip --out-dir .

# Assembly 数据库
$ esearch -db assembly -query 'txid670[Organism] AND latest[filter] NOT anomalous[filter]' | efetch -format xml > vpa.assembly.xml
# Assembly 数据库对应的 Biosample 去除冗余的样本信息
$ for i in $(xtract -input vpa.assembly.xml -pattern DocumentSummary -element AssemblyAccession BioSampleAccn | awk '{print $2}' | awk '!a[$0]++'); \
> do efetch -db Biosample -id $i > biosample/assembly/$i.xml; sleep 5; \
> done
```
