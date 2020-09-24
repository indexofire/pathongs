# Nanopore 基因组测序

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! "Abstract"


## 案例1.病毒比对获得基因组

针对一个1.5k～1.8k左右的Dengue病毒扩增子测序文库进行的 Nanopore 测序数据处理。共对24个样本进行文库构建测序，使用NBD104,114作为barcodes。

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

### 1.测序结果扫描

非必须，初步了解目标扩增子在测序结果中的分布，如果有大量人源reads，可能扩增效率不佳，需要改换引物重新扩增。

```bash
# cenrifuge 对 reads 进行扫描
(centrifuge)$ export DB=$HOME/dbs/centrifuge/b_compressed+h+v
(centrifuge)$ ls 00_raw/*.fastq | parallel --max-args=1 centrifuge -x $DB -U {1} -S 01_centrifuge/{1}.result --report-file 01+centrifuge/{1}.reprot
```

### 2.去除接头

使用 guppy 做 demultiplex 和 basecalling 时没有做trim barcodes，需要使用 porechop 等软件去除 barcodes 序列。

```bash
(nanopore)$ ls 00_raw/*.fastq | parallel --max-args=1 porechop -i {1} -o 02_porechop/{1}
```

### 3.去除低质量序列和短片段

对于 Nanopore 生成的碱基质量做一便过滤，这里使用reads Q值为9作为cut-off，reads 长度超过1200的保留。

```bash
(nanopore)$ ls 02_porechop/*.fastq | parallel --max-args=1 NanoFilt -q 9 -l 1200 {1} | gzip > 03_nanofilt/trimed.fastq.gz
```

### 4.去除引物区域

去除5'和3'端的引物区域。如果reads比对时起始比对的位置落在引物扩增区，则去除。

```bash

```


## 2.细菌组装获得基因组

案例使用霍乱弧菌gDNA做为起始样本，进行de novo测序。获得 Nanopore reads 后直接使用组装软件对nanopore数据进行组装获得基因组。

常用的组装软件：

- flye
- canu
- shata

**使用 flye**

```bash
# flye 进行基因组组装 nanopore 测序的霍乱基因组 reads
(nanopore)$ flye --nano-raw trimed.fq.gz -g 4m -t 160 -o 04_flye
(nanopore)$ ls 04_flye
assembly.fasta    # 组装结果，含contigs的fasta格式文件
```

**使用 canu**

```bash
(nanopore)$ canu
```

**使用 miniasm**

```bash
(nanopore)$
```

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

### 构建流程

**采用shell脚本**

```bash
#!/bin/sh
# pipeline

# 过滤获得组装reads数据
for i in barcode{01..24}
  do mkdir ${i}/02_porechop ${i}/03_nanofilt
  ls ${i}/*.fastq | parallel --max-args=1 porechop -i {1} -o {//}/02_porechop/{1}
  ls ${i}/02_porechop/*.fastq | parallel --max-args=1 NanoFilt -q 9 {1} | gzip > ${i}/03_nanofilt/${i}_trimed.fq.gz
done

# 组装成contigs
ls barcode*/03_nanofilt/*.fq.gz | parallel --max-args=1 flye --nano-raw {1} -g 4m -o {//}/../04_flye
```

**采用snakemake**

```bash

```

## Reference




```bash
#!/bin/sh

DB=$HOME/dbs/centrifuge/b_compressed+h+v
CPU=160

for i in barcode{01..24}; do
    mkdir ${i}/00_raw ${i}/01_centrifuge ${i}/02_porechop ${i}/03_nanofilt ${i}/04_mapping
    mv ${i}/*.fastq ${i}/00_raw
    ls ${i}/00_raw/*.fastq | parallel --max-args=1 centrifuge -x $DB -U {1} -S ${i}/01_centrifuge/{1}.result --report-file ${i}/01_centrifuge/{1}.reprot
    cd ${i}/00_raw
    ls *.fastq | parallel --max-args=1 porechop -i {1} -o ../02_porechop/{1}
    cd ../02_porechop
    ls *.fastq | parallel --max-args=1 NanoFilt -q 9 -l 1200 {1} | gzip > ../03_nanofilt/${i}_trimed.fastq.gz
    cd ../../
    bwa index ref.fa
    bwa mem -t $CPU -R $(echo "@RG\tID:${i}\tSM:${i}") ref.fa ${i}/03_nanofilt/${i}_trimed.fastq.gz | samtools view -bS - | samtools sort -@$CPU -O bam -o ${i}/04_mapping/${i}_sorted.fastq.gz
    bwa index ${i}/04_mapping/${i}_sorted.fastq.gz
done
```


```bash
mkdir $1/00_raw $1/01_centrifuge $1/02_porechop $1/03_nanofilt $1/04_mapping
mv $1/*.fastq $1/00_raw
centrifuge -x $DB -U $ -S ${i}/01_centrifuge/{1}.result --report-file ${i}/01_centrifuge/{1}.reprot
cd ${i}/00_raw
ls *.fastq | parallel --max-args=1 porechop -i {1} -o ../02_porechop/{1}
ls *.fastq | parallel --max-args=1 NanoFilt -q 9 -l 1200 {1} | gzip > ../03_nanofilt/${i}_trimed.fastq.gz
cd ../../
bwa index ref.fa
bwa mem -t $CPU -R $(echo "@RG\tID:${i}\tSM:${i}") ref.fa ${i}/03_nanofilt/${i}_trimed.fastq.gz | samtools view -bS - | samtools sort -@$CPU -O bam -o ${i}/04_mapping/${i}_sorted.fastq.gz
bwa index ${i}/04_mapping/${i}_sorted.fastq.gz



```



```bash
#!/bin/sh
CPU=160
HEADER=$(zcat $1 | head -n 1)
ID=$(echo $HEADER | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
SM=$(echo $HEADER | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

bwa mem -M -t $CPU -v 3 -R $(echo "@RG\tID:$ID\tSM:$ID"_"$SM\tLB:$ID"_"$SM\tPL:ILLUMINA") \
ref.fa $1 $2 | samtools view -bS - | samtools sort -@$CPU -O bam -o $ID_$SM_mapped.bam
```



```bash
# snakemake


```


## Nanopore 数据直接组装的应用

虽然目前如 flye 等组装软件对于细菌基因组的组装准确度已经很高了，但是由于 Nanopore 本身技术特点原因，还是不可能像2代测序一样获得精准的基因组序列的。

- 细菌纯 Nanopore 数据进行组装不能胜任：基于SNP分析，比如MLST，MLVA，基于SNP构建进化树等。要进行这些分析需要结合二代数据。
- 可以做的：基于基因和基因组结构的分析，比如细菌基因组结构共线性分析，毒力基因、耐药基因的扫描等

## Reference

- https://timkahlke.github.io/LongRead_tutorials
