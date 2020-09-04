# 下载公共数据库中基因组数据

---

@indexofire update: 2020.4.21

```bash
# 构建环境
$ conda create -n dlngs
$ conda activate dlngs
```

## 下载 Assembly 数据

```bash
(dlngs)$ conda install entrez-direct
```


## 从 NCBI 下载 sra 数据

首先要安装 sra-tools 工具

```bash
(dlngs)$ conda install sra-tools
```

- prefetch 下载
- fastq-dump
- fasterq-dump
- vdb-config 设置 sra-tools 工具的参数

::: Note

目前 NCBI 已将 sra 数据放到 AWS 和 GCP 云服务器端，过去用 aspera 快速下载的方式已经无法使用，只能用 https 模式下载。如果你在国外，或者有快速的代理服务器，还可以获得比较优秀的下载速度。

从 AWS 或 GCP 下载。国内访问 GCP 需要代理，而AWS的S3服务器也不定时抽风。


## 从 EBI 下载 FastQ 数据

EBI ENA 支持直接下载 fastq 格式的基因组数据，可以用 aspera 加速下载。

```bash
(dlngs)$ conda install aspera-cli
```

**使用 ascp 下载 fastq 数据示例**

```bash
(dlngs)$ ascp -QT -l 100m -P33001 -i $CONDA_PREFIX/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR100/000/SRR1001720/SRR1001720_1.fastq.gz .
```

**使用 ftp 模式下载**

```bash
(dlngs)$ wget ftp://ftp.ebi.ac.uk/vol1/fastq/SRR100/000/SRR1001720/SRR1001720_1.fastq.gz
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
