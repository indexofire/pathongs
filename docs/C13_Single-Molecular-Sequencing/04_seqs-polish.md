# 数据清洗

---

!!! Abstract "内容介绍"
    因为Nanopore序列准确性的特点，需要进行数据清洗已完成后续的分析工作，本节以 Minion 为例介绍纳米孔测序下机数据的清洗操作。

## 1. porechop

去除接头

```bash
$ porechop -i minion_brown_metagenome/brown_metagenome.template.fasta > brown_metagenome.template.chopped.fasta

```
