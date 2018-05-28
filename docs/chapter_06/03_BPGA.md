# 3.2 [BPGA][]

[BPGA][] 是由印度 [CSIR-Indian Institute of Chemical Biology](http://www.iicb.res.in/) 研究人员开发的用于细菌 Pangenome 分析工具流。

[BPGA][] 依赖第三方软件，用 [USEARCH][] 构建聚类，用 [Muscle][] 进行多重序列比对。用 [gnuplot][] 和 [ghostscript][] 生成结果图片或pdf文件。

## 使用

**输入数据**

* .faa 蛋白质序列文件
* .pep.fsa HMP蛋白数据文件
* .gbk/gb [Genbank][]蛋白质序列文件

**使用方法**

BPGA 采用了类似phylip式的交互式运行方式。

```
$ BPGA-Version-1.3
```

运行程序后，会进入菜单选择界面，一般可以分2步走：

1. 按 1 选择输入文件格式
2. 按 2 进入 Pangenome 分析

按 0 退回上一级菜单

**结果文件**

生成的数据包括根目录下的文件和Results、Sequences 和 Supporting_files 3个子文件夹内结果文件

Results文件夹内主要包含分析结果生成的图：

| 结果文件 | 数据内容 |
| -------- | -------- |
| COG_DISTRIBUTION_DETAILS.pdf | core, accessory, unique 3种类型基因在COGS各个分类中的分布比值 |
| COG_DISTRIBUTION.pdf | core, accessory, unique 3种类型基因在COGS中总的分布比值 |
| Core_Pan_Dot_Plot.pdf | pangenome 和 coregenome 点阵趋势图 |
| Core_Pan_Plot.pdf | pangenome 和 coregenome boxplot趋势图 |
| curve.xls | pangenome 和 coregenome 基因数计算公式和物种的基因估计值，文本格式 |
| Default_Core_Pan_Plot.pdf | pangenome 和 coregenome 趋势图 |
| Histogram.pdf | 各个分析菌株所含基因家族的直方图分布 |
| KEGG_DISTRIBUTION_DETAILS.pdf | pdf 格式的 KEGG 详细分类基因的 core, accessory 和 uniq genes 数量分布 |
| KEGG_DISTRIBUTION.pdf | pdf 格式的 KEGG 主要分类基因的 core, accessory 和 uniq genes 数量分布|
| New_Genes_Plot.pdf | 分析菌株所含 New genes 数量分布 |
| Pan_phylogeny.pdf | pdf 格式的分析菌株提供发生树图|
| Pan_phylogeny.svg | svg 格式的分析菌株系统发生树图 |
| stats.xls | 所有分析菌株含有的 core, accessory, unique 和 exclusively absent 基因。虽然后缀为.xls，文本格式 |

Sequences 文件夹内包含：

| 结果文件 | 数据内容 |
| -------- | -------- |
| accessory_seq.txt | accessory 基因的序列 |
| core_genes_with_atypical_GC_content.txt | 进行了 atypical GC 分析后产生， core 基因的序列 |
| core_seq.txt | uniq 基因的序列 |
| exclusively_absent_seq.txt | exclusively absent 基因的序列 |
| REPSEQ_ACCESSORY.txt | |
| REPSEQ_CORE.txt | |
| REPSEQ_UNIQUE.txt | |
| unique_genes_with_atypical_GC_content.txt | 进行了 atypical GC 分析后产生， uniq 基因的序列 |
| unique_seq.txt | uniq 基因的序列 |

Supporting_files 文件夹包含：

| 结果文件 | 数据内容 |
| -------- | -------- |
| ACCESSORY_COG_hits3.txt | |
| ACCESSORY_kegg_hits3.txt | |
| Cog_Category1.txt | |
| CORE_COG_hits3.txt | |
| core_default.txt | |
| core_genome.txt | |
| CORE_kegg_hits3.txt | |
| DATASET.xls | |
| histogram.txt | |
| kegg_accessory_id.txt | |
| kegg_accessory_out.txt | |
| kegg_core_id.txt | |
| kegg_core_out.txt | |
| Kegg_count_details1.txt | |
| kegg_histogram1.txt | |
| kegg_unique_id.txt | |
| kegg_unique_out.txt | |
| list | |
| Major_Cog_Category1.txt | |
| matrix.txt | |
| new_genes_count.txt | |
| pan_default.txt | |
| pan_genome.txt | |
| PAN_PHYLOGENY_MOD.nwk | |
| PAN_PHYLOGENY_MOD.ph | |
| PAN_PHYLOGENY.ph | |
| plots_default.plt | |
| plots.plt | |
| u_clusters.txt | |
| UNIQUE_COG_hits3.txt | |
| UNIQUE_kegg_hits3.txt | |
