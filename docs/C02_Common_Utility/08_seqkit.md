# Seqkit 好用的命令行序列处理工具套件

Seqkit 是由国人开发的序列处理工具

!!! note "例子1"
    对注释的多拷贝基因提取序列后进行序列比对。

```bash
# 用prokka注释基因组拼接序列
$ prokka --rfam --prefix sample contig.fa

# 用 seqkit 提取多拷贝基因 ncRNA Qrr 的序列
$ grep 'Qrr' sample.gff | for i in $(awk '{print $9}'); do seqkit grep -p ${i:3:14} sample.ffn >> qrr.fa; done

# mafft 序列比对
$ mafft 1.fa > 1.mafft

# raxml 构建进化树
$ raxmlHPC -f a -p 12345 -x 12345 -# 1000 -m GTRGAMMA -s 1.mafft -n 1
```
