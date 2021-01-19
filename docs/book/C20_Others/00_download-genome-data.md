# 命令行下载公共数据库中基因组数据

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "本节简介"
    简介使用命令行工具从公共数据库中下载各种基因组数据。

共公共数据库中下载微生物基因组数据，最常用的从 NCBI 的 Assembly, Genbank 等数据库下载 fasta/genbank 格式的数据，或者从 SRA 下载测序原始数据。

!!! Warning "SRA数据下载方式更新"
    目前 NCBI 已将新提交的 SRA 数据放到 AWS 和 GCP 云服务器端，NCBI服务器只保留了部分老的 SRA 数据。过去用 aspera 快速下载新数据的方式已经无法使用，如果调用 prefetch 只能用 https 模式下载，国内速度会比较慢。如果你在国外，或者有快速的代理服务器，还可以获得比较优秀的下载速度。如果安装了 AWS 和 GCP 的命令行工具下载，国内访问 GCP 需要代理，而 AWS 的 S3 服务器也不定时抽风，基本上也要靠代理。

首先创建一个虚拟环境，安装所需的下载工具。

```bash
# 构建环境
$ conda create -n dlngs
$ conda activate dlngs
```

## 1.下载 Assembly 数据

```bash
# 下载 NCBI entrez-direct 工具
(dlngs)$ conda install entrez-direct
#
(dlngs)$ esearch -query | efetch
```

## 2.下载 SRA 数据

首先要安装 sra-tools 工具

```bash
(dlngs)$ conda install sra-tools
```

- prefetch 下载
- fastq-dump
- fasterq-dump
- vdb-config 设置 sra-tools 工具的参数

## 从 EBI 下载 FastQ 数据

EBI ENA 支持直接下载 fastq 格式的基因组数据，可以用 aspera 加速下载。

```bash
(dlngs)$ conda install -c hcc aspera-cli
```

**使用 ascp 下载 fastq 数据示例**

```bash
(dlngs)$ ascp -QT -l 100m -P33001 -i $CONDA_PREFIX/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR100/000/SRR1001720/SRR1001720_1.fastq.gz .
```

**使用 ftp 模式下载**

```bash
(dlngs)$ aria2c ftp://ftp.ebi.ac.uk/vol1/fastq/SRR100/000/SRR1001720/SRR1001720_1.fastq.gz
```

如果我们希望批量下载，可以用下面的脚本。现将要下载的SRA数据的Accession保存成list.txt文件，然后用`ena-dl.sh list.txt output`下载到output目录

```bash
#!/bin/sh
# ena-dl.sh list.txt

for i in $(cat $1);
do
  ascp -QT -l 100m -P33001 -i $CONDA_PREFIX/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/${i:0:6}/00${i:0-1:1}/$i/${i}_1.fastq.gz $2/${i}_1.fastq.gz;
  ascp -QT -l 100m -P33001 -i $CONDA_PREFIX/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/${i:0:6}/00${i:0-1:1}/$i/${i}_2.fastq.gz $2/${i}_2.fastq.gz;
done
```

如果是SE测序，那么只有_1.fastq.gz文件可以被下载。如果是6位长度的SRA文件，则不需要000的末尾编号。
