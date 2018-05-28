# get_homologues

[get_homologues](http://eead-csic-compbio.github.io/get_homologues) 是一个用 perl 语言编写的用来鉴定细菌 core-genomes 和 pan-genomes 的开源工具。

## 分析案例

我们以 Bacillus cereus 基因组 assembly 数据库为例分析该物种的 Pangenomics。

```bash
# 抓取 gbk 格式的 assembly 文件在 NCBI FTP 上的下载路径
$ esearch -db assembly -query "Bacillus cereus[ORGN] AND latest[SB]" | \
> efetch -format docsum | \
> xtract -pattern DocumentSummary -element FtpPath_RefSeq | \
> awk -F"/" '{print $0"/"$NF"_genomic.gbff.gz"}' > bcereus.path

# 用 wget 工具下载
$ wget --limit-rate 300k --no-passive-ftp -i bcereus.path

# get_homologues 分析
$ get_homologues.pl -d . -n 40
```

get_homologues 可以默认使用 blast 进行序列相似性搜索，为了加快速度，也可以使用 diamond。

```bash
# -X 调用 diamond 进行序列相似性搜索
$ get_homologues.pl -d gbk_folder -n 40 -X
```






## Reference

1. Bacterial Pangenomics, Methods and Protocols, Chapter14
2. [get_homologues manual](http://eead-csic-compbio.github.io/get_homologues/manual/manual.html)

