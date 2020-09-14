# Abricate

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 安装

```bash
$ conda create -n abricate
$ conda activate abricate
(abricate)$ conda install abricate
```

## 使用

```bash
# 查看支持的数据库，默认采用resfinder
(abricate)$ abricate --list
# 扫描 fasta
(abricate)$ abricate genome.fasta

# 选择其他数据库 vfdb
(abricate)$ abricate -db vcdb genome.fasta

# 批量扫描生成结果
(abricate)$ for i in *.fa; do abricate --db card --quiet $i > ${i%.fa}.result; done
(abricate)$ abricate --summary *.result > result
(abricate)$ cat result
```


1. Abricate: https://github.com/tseemann/abricate
