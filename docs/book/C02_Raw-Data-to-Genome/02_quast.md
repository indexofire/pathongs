# Quast 质量分析

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 安装

```bash
# Quast 依赖 python2.x，因此新建一个 virtualenv
$ conda create -n quast
$ conda activate quast
(quast)$ conda install quast
```

## 使用

```bash
# -o 制定输出目录 contig.fa 为拼接结果序列
$ quast -o outdir assembly/contig.fa
# 浏览器打开报告
$ chromium outdir/report.html
```

- N50: Reads拼接后会获得一些不同长度的Contigs。将所有的Contig长度相加，能获得一个Contig总长度。然后将所有的Contigs按照从长到短进行排序，将Contig按照这个顺序依次相加，当相加的长度达到Contig总长度的一半时，最后一个加上的Contig长度即为Contig N50。举例：Contig 1+Contig 2+ Contig 3 +Contig 4=Contig总长度*1/2时，Contig 4的长度即为Contig N50。Contig N50可以作为基因组拼接的结果好坏的一个判断标准。
- N75: 与N50方法一致，不过长度为超过Contigs总长度75%时的contig序列长度。
- L50: 为获得N50时contigs数量。
- L75: 为获得N75时contigs数量。
