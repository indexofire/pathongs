# centrifuge

---

!!! Abstract "内容简介"
    Centrifuge 是一个序列taxnomic分类软件，可以用来做宏基因组分析，也常用来做微生物测序序列污染鉴定。

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

对测序数据或者fasta序列进行物种鉴定扫描。

```bash
# 扫描 fasta 序列
(centrifuge)$ centrifuge -x $HOME/dbs/centrifuge/p_compressed+h+v -U example.fa --report-file report.txt -S results.txt
# 扫描 fastq 序列
(centrifuge)$ centrifuge -x $HOME/dbs/centrifuge/p_compressed+h+v -U S1_R1.fastq --report-file S1_R1-report.txt -S S1_R1-results.txt
# 扫描多个 fastq 序列
(centrifuge)$ for i in *.fq.gz; do zcat $i >> total.fq.gz; done
(centrifuge)$ centrifuge -x $HOME/dbs/centrifuge/p_compressed+h+v -U total.fq.gz
```

对序列结果进行分析。可以用 --report-file 参数生成的报告文件，里面将各个物种进行了汇总。或者用`centrifuge-kreport`工具将reads扫描的结果默认文件 centrifuge_report.tsv 转换成 kraken 格式的物种分类树状结果。这个数据结果可以用[pavian](https://github.com/fbreitwieser/pavian)等进行数据可视化。

```bash
(centrifuge)$ centrifuge-kreport -x $HOME/dbs/centrifuge/p_compressed+h+v centrifuge_report.tsv > kraken.result
```
