# MLST 与基因组数据

---

!!! Abstract "内容简介"


| tools | fasta | fastq | gz 压缩 |
| ----- | ----- | ----- | ------- |
| mlst | YES | NO | |
| srst2 | | YES | |


### mlst

**安装**

```bash
# conda 安装
$ conda create -n mlst
$ conda activate mlst
(mlst)$ conda install mlst
```

**更新数据库**

```bash
# 更新数据库
(mlst)$ cd $CONDA_PREFIX/scripts
(mlst)$ ./mlst-download_pub_mlst | bash
(mlst)$ mv ../db/pubmlst ../db/pubmlst.old
(mlst)$ cp -R pubmlst ../db
# 建立新数据库索引
(mlst)$ ./mlst-make_blast_db
```

**自定义数据库**

```bash


```

**检索fasta数据的MLST型别**

```bash

```

### srst2

对文件名命名格式有要求，miseq下机数据的格式，或者_1.fastq.gz, \_2.fastq.gz格式。

```bash
(srst2)$ srst2 --input_pe
```
