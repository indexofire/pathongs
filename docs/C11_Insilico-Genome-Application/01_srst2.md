# SRST2

---

[SRST2](https://github.com/katholt/srst2) 工具可以对 reads 进行扫描，分析其MLST型別以及基因部分情况。

## 1. 安装

```bash
# 使用 conda 安装 SRST2 安装包
$ conda create -n srst2 srst2
$ conda activate srst2
```

## 2. 使用

需要注意的是fastq文件名称

```bash
# 下载物种 MLST 数据库
$ getmlst.py --species "Salmonella enterica"

# 根据测序reads扫描沙门菌MLST型别
$ srst2 --input_pe S1_1.fastq.gz S1_2.fastq.gz --output S1 \
> --log --mlst_db Salmonella_enterica.fasta \
> --mlst_definitions senterica.txt --mlst_delimiter _

# 根据测序reads获得耐药基因
$ srst2 --input_pe
```

## 3. 实例

下载 Bacillus anthracis SRA 数据库miseq测序平台的基因组测序数据。

```bash
# edirect 工具生成 sra 的 runinfo 信息列表
# 过滤了测序数据过小和平均读长过小的数据
$ esearch -query '"Bacillus anthracis"[Organism] AND \
> "miseq"[All Fields] AND ("biomol dna"[Properties] AND \
> "strategy wgs"[Properties])' -db sra | efetch -format runinfo \
> -db sra | awk '/^[SDE]RR/' | \
> awk -F',' '{if($8>150  && $16=="PAIRED" && $20 =="Illumina MiSeq") print $1}' | \
> prefetch -v
```

将 SRA 格式的数据转换成 fastq 格式。

```bash
# 批量转化成 fastq.gz 格式文件
$ parallel "fastq-dump --split-files --gzip --outdir fastq" ::: *.sra
# 确认数据均是 paired end 测序，如果有 single end，将其分离到不同的目录中
$ ls -l *.fastq.gz | awk -F'_' '{print $1}' | awk '{print $9}' | uniq -u
# 输出的编号即为非 paired end 测序数据。

# bioawk 统计GC含量分布
$ parallel "bioawk -c fastx 'BEGIN{n=0;q=0}{n+=gc(\$seq);q+=meanqual(\$seq)}END{print \$name,n/NR,q/NR}' \
> >> gc_result.txt" ::: *.fastq.gz
```

srst2 对数据进行扫描，获得MLST型别，毒力基因

```bash
# 准备目录结构
$ mkdir -p mlst vfdb ardb result

# 准备 MLST 数据库，将下载的文件放入 mlst 目录内
# 由于 B. anthracis 没有公共 MLST 数据库，这里用近源的 B. cereus MLST 数据库
# 数据来分析 B. anthracis
$ getmlst.py --species "Bacillus cereus"

# 准备毒力基因数据文件，将文件放入 vfdb 目录中
$ cd vfdb
$ wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
$ gunzip VFDB_setB_nt.fas.gz
$ python path/srst2/database_clustering/VFDBgenus.py --inflie VFDB_setB_nt.fas --genus Bacillus
$ cd-hit -i Bacillus.fsa -o Bacillus_cdhit90 -c 0.9 > Bacillus_cdhit90.stdout
$ python path/srst2/database_clustering/VFDB_cdhit_to_csv.py --cluster_file Bacillus_cdhit90.clstr \
> --infile Bacillus.fsa --outfile Bacillus_cdhit90.csv
$ python path/srst2/database_clustering/csv_to_gene_db.py -t Bacillus_cdhit90.csv \
> -o Bacillus_VF_clustered.fasta -s 5

# 批量处理
$ for i in $(ls *.fastq.gz | awk -F'_' '{print $1}' | uniq); do srst2 \
> --input_pe ${i}_1.fastq.gz ${i}_2.fastq.gz \
> --mlst_db ../mlst/Bacillus_cereus.fasta \
> --mlst_definitions ../mlst/bcereus.txt \
> --mlst_delimiter _ \
> --gene_db Bacillus_VF_clustered.fasta \
> --log --threads 40 --output ../result; \
> done

# 耐药基因分析
$ for i in $(ls *.fastq.gz | awk -F'_' '{print $1}' | uniq); do srst2 \
> --input_pe ${i}_1.fastq.gz ${i}_2.fastq.gz \
> --log --output ../result \
> --gene_db path/srst2/datab/ARGannot.r1.fasta; \
> done
```
