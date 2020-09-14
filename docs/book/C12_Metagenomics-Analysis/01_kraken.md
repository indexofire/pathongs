# Kraken

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    Kraken 是一个序列taxnomic分类软件，可以用来做宏基因组分析，也常用来做微生物测序序列污染鉴定。

## 1. 安装

[kraken2](https://www.ccb.jhu.edu/software/kraken2/) 是目前最新的版本。

```bash
# 建立conda环境
$ conda create -n kraken kraken2 bracken=2.2
$ conda activate kraken
(kraken)$ wget https://www.ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz
# 下载比对数据库
(kraken)$ tar zxvf minikraken2_v2_8GB.tgz -C $HOME/dbs/minikraken2
```

## 2. 使用

使用kraken2需要用`--db`指向其数据库，比如微生物鉴定常用minikraken2数据库。用`--report`参数生成结果报告。

```bash
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report.txt assembly.fna > result.kraken
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report2.txt --fastq-input --paired S1_R1.fq.gz S1_R2.fq.gz > result2.kraken
(kraken)$ wget https://ccb.jhu.edu/software/bracken/dl/minikraken2_v2/database200mers.kmer_distrib
(kraken)$ mv database200mers.kmer_distrib $HOME/dbs/minikraken2
(kraken)$ bracken -d $HOME/dbs/minikraken2 -i kraken2.report -o bracken.species.txt -l S
```

## 3. 示例
