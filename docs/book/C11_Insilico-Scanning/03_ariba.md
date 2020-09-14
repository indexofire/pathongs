# ARIBA

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

[Ariba](https://github.com/sanger-pathogens/ariba) 工具是由 Sanger 测序中心病原微生物开发组开发的用于Assembly数据扫描的工具。主要用于耐药基因的扫描，也可以用来对特定基因扫描（如MLST管家基因）。

## 1. 安装

```bash
# bioconda 安装包版本较最新版落后较多，通过 pip 从pypi获取最新版。通过 pip 安装就要手动安装依赖软件，这里还是以 conda 安装为例。
$ conda create -n ariba python=3.6
(ariba)$ conda install bowtie2 cd-hit mummer
(ariba)$ pip install ariba
```

## 2. 使用
