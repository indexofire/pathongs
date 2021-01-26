# Abricate

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容介绍"
    本节介绍一个非常好用的细菌基因扫描工具: abricate，这是由澳大利亚著名生物信息学家[@Torsten Seemann](https://twitter.com/torstenseemann) 开发的细菌毒力基因和耐药基因扫描工具。

abricate 是一个很方便的工具，可以对细菌基因组 fasta 数据进行 blast 扫描，获得其基因组携带那些毒力基因和耐药基因结果。软件包含数据库，不需要自己再额外下载构建，数据库包括：

- vfdb: 基于vfdb毒力基因数据库
- argannot: 基于argannot耐药数据库
- plasmidfinder: 基于plasmidfinder用于质粒扫描
- ecoli_vf: 基于大肠杆菌毒力基因
- megares: 基于megares耐药数据库
- resfinder: 基于resfinder耐药数据库
- ecoh: 用于大肠杆菌O/H抗原分析
- ncbi: 基于NCBI的耐药数据库
- card: 基于加拿大CARD耐药数据库

## 安装

使用版本：

- abricate: 1.0.1

通过 conda 安装 abricate，无需安装其他依赖。

```bash
# 新建一个虚拟环境运行 abricate
$ conda create -n abricate
$ conda activate abricate
# 安装 abricate
(abricate)$ conda install abricate
```

## 使用

abricate 内建多个数据库，可以用`--list`查看，要使用哪个数据库，则采用`-db`参数加上具体数据库名称即可。

```bash
# 查看支持的数据库，默认采用resfinder
(abricate)$ abricate --list
DATABASE        SEQUENCES       DBTYPE  DATE
plasmidfinder   460     nucl    2020-Apr-19
vfdb    2597    nucl    2020-Apr-19
argannot        2223    nucl    2020-Apr-19
ecoli_vf        2701    nucl    2020-Apr-19
megares 6635    nucl    2020-Apr-19
resfinder       3077    nucl    2020-Apr-19
ecoh    597     nucl    2020-Apr-19
ncbi    5386    nucl    2020-Apr-19
card    2631    nucl    2020-Apr-19

# 扫描 fasta
(abricate)$ abricate genome.fasta

# 选择其他数据库 vfdb
(abricate)$ abricate -db vfdb genome.fasta

# 批量扫描生成结果
(abricate)$ for i in *.fa; do abricate --db card --quiet $i > ${i%.fa}.result; done
(abricate)$ abricate --summary *.result > result
(abricate)$ cat result

# 用parallel 加速并行结果
(abricate)$ ls *.fa | parallel --max-args=1 abricate --db vfdb -q {1} > vfdb.result
(baricate)$ less vfdb.result
```

!!! Warning "注意事项"
    abricate 过滤条件中对 DNA identity 和 coverage 默认设置都是 80%。

## 自定义

可以用abricate扫描自行构建的基因数据库。

## Reference

1. Abricate: https://github.com/tseemann/abricate
