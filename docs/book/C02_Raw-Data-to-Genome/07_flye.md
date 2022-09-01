# Nanopore测序数据组装

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    本节介绍我们实验室常用的对细菌纳米孔测序数据进行直接组装的方法。随着Guppy引入SUP Basecalling算法，获得序列一致性达到～95%的水平，针对该质量水平的raw data进行直接组装获得细菌基因组已经可以获得一个非常优秀的组装结果，可以进行耐药基因分析，甚至可以获得准确的MLST型别。虽然进化分析的影响还需要更进一步的评估才能确定，但相信随着接下来Oxford Nanopore公司推出的Q20试剂，以及Q30实现方法，将会对细菌基因组测序带来更为准确和便捷的模式。我们实验室进行测试后将会尽快更新细菌组装方案。

## Flye + racon + medaka 流程

**使用软件:**

- Flye v2.9
- racon v1.4.0
- medaka v1.5.0
- filtlong
- pyco
- mummer4

**基本分析流程**

### 安装环境

```bash
$ conda create -n ont-assembly
$ conda activate ont-assembly
(ont-assembly)$ mamba install medaka
(ont-assembly)$ mamba install flye racon filtlong
(ont-assembly)$ mamba install mummer4
```

### 分析流程

**数据预处理**

本例子采用的数据是在Gridion的R9.4芯片上测序下机的细菌基因组DNA测序数据。样本基因组gDNA没有进行片段化处理，文库构建采用SQK-RBK004进行。Basecalling 采用Guppy v5.0.17 SUP Model。对获得的数据首先进行QC质检：

```bash
# 本教程的数据
(ont-assembly)$ ls data/
FAQ89450_pass_barcode01_f72491a8_0.fastq
FAQ89450_pass_barcode01_f72491a8_1.fastq
FAQ89450_pass_barcode01_f72491a8_2.fastq
...

(ont-assembly)$
(ont-assembly)$ filtlong
(ont-assembly)$ flye --nano-hq data -t 16 -g 5m -o flye-output
--nano-hq  表示所用的纳米孔测序数据错误率为小于5%的数据
-g 5m      表示组装目标物种的基因组大小约5M
-t 16      表示采用16个线程进行基因组组装
```

## shasta + pepper + margin


## wt


## 常用流程分析结果比较

-
-
-

我们对本例使用的数据集进行了结果比较
