# 混合组装

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

结合 illumina 测序和 PacBio/Nanopore 测序数据，可以较容易获得细菌基因组的完成图。常见的混合组装软件：

- spades
- canu


## 方法

### 1. spades 组装

```bash
(denovo)$ spades.py
```

### 2. unicycler 组装

```bash
$ conda activate nanopore
(nanopore)$ unicycle -1 -2 -l
```
