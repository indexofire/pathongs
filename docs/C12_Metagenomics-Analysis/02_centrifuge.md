# centrifuge

---

!!! Abstract "内容简介"
    Centrifuge 是一个序列taxnomic分类软件，可以用来做宏基因组分析，也可以用来做微生物测序序列污染鉴定。

- 版本: v1.0.3 beta

## 安装

```bash
# 新建一个虚拟环境，安装 centrifuge
$ conda create -n centrifuge centrifuge
# 进入虚拟环境
$ conda activate centrifuge
# 下载数据库
(centrifuge)$ aria2c ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz
# 解压缩数据库文件
(centrifuge)$ tar zxvf p_compressed+h+v.tar.gz -C $HOME/dbs/centrifuge
```

## 使用

```bash
# 扫描 fasta 序列
(centrifuge)$ centrifuge -x $HOME/dbs/centrifuge/p_compressed+h+v -U example.fa --report-file report.txt -S results.txt

# 扫描 fastq 序列
(centrifuge)$ centrifuge -x $HOME/dbs/centrifuge/p_compressed+h+v -U S1_R1.fastq --report-file S1_R1-report.txt -S S1_R1-results.txt
```
