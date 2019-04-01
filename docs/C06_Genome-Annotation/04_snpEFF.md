# snpEFF 使用

snpEFF 是许多 variant calling pipeline 中常用的 vcf 注释工具。

## 安装

```bash
$ conda install snpEFF

```

```bash
# 将 Listeria monocytogenes 涉及的数据库下载
$ snpEFF databases | awk '/^Listeria_monocytogenes'/ {print $4}' | awk '!a[$1]++{print}' | xargs wget
```
