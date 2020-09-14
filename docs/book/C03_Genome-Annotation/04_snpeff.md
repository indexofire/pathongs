# snpEFF 使用

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

snpEFF 是许多 variant calling pipeline 中常用的 vcf 注释工具。

## 安装

```bash
$ conda install snpEFF
```

## 使用

```bash
# 将 Listeria monocytogenes 涉及的数据库下载
$ snpEFF databases | awk '/^Listeria_monocytogenes'/ {print $4}' | awk '!a[$1]++{print}' | xargs wget
```
