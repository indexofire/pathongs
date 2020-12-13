# Kraken 使用简介

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    Kraken 是一个序列taxnomic分类软件，可以用来做宏基因组分析，也常用来做微生物测序序列污染鉴定。

本文使用软件版本：

- [kraken2](https://www.ccb.jhu.edu/software/kraken2/): v2.1.0
- [braken](https://ccb.jhu.edu/software/bracken/) : v2.5.0

## 1. 安装

[kraken2](https://www.ccb.jhu.edu/software/kraken2/) 是目前最新的版本，这里使用 conda 方式安装。

```bash
# 新建 conda 环境，安装 kraken2 和 bracken
$ conda create -n kraken kraken2 bracken
$ conda activate kraken
# 下载 minikraken 数据库
(kraken)$ wget https://www.ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz
# 将数据库解压缩至指定目录，这里放在$HOME/dbs/minikraken2
(kraken)$ tar -zxf minikraken2_v2_8GB.tgz -C $HOME/dbs/minikraken2
```

## 2. 使用

使用kraken2需要用`--db`指向其数据库，比如微生物鉴定常用 minikraken2 数据库。用`--report`参数生成结果报告。

```bash
# 对宏基因组组装结果 assembly.fna 序列进行扫描
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report.txt assembly.fna > result.kraken
# 对测序reads数据进行扫描
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report2.txt --fastq-input --paired S1_R1.fq.gz S1_R2.fq.gz > result2.kraken
```

bracken 是采用贝叶斯算法的 kraken 丰度评估软件。

```bash
# 采用kmer方式对已建立kmer索引的物种数据库进行比对，获得序列对应物种
(kraken)$ wget https://ccb.jhu.edu/software/bracken/dl/minikraken2_v2/database200mers.kmer_distrib
(kraken)$ mv database200mers.kmer_distrib $HOME/dbs/minikraken2
(kraken)$ bracken -d $HOME/dbs/minikraken2 -i kraken2.report -o bracken.species.txt -l S
```
