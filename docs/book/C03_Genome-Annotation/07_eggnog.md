# EggNOG

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

- [EggNOG](http://eggnogdb.embl.de/#/app/home)
- [EggNOG-mapper](https://github.com/jhcepas/eggnog-mapper)

## 本地安装 EggNOG-mapper

```bash
# 安装
$ conda create -n eggnog
$ conda activate eggnog
(eggnog)$ conda install eggnog-mapper
# 我们只关心细菌和病毒的基因组注释，因此只下载这2个数据库
# 各个单独数据库参见 http://beta-eggnogdb.embl.de/download/eggnog_4.5/hmmdb_levels
(eggnog)$ cd PATH/TO/DATA
(eggnog)$ download_eggnog_data.py -yq --data_dir `pwd` bact
```

## 使用

emapper 的默认算法是hmm，如果想用diamond来做cluster，需要添加`-m diamond`参数

```bash
# HMM 注释
# 针对细菌基因组数据库，使用40个CPU内核，数据库运行在内存中
(eggnog)$ emapper.py -i mygenome.fa -d bact --data_dir PATH/TO/DATA \
> --cpus 40 --usemem -o mygenome_bact_result
# 针对病毒基因组数据库
(eggnog)$ emapper.py -i mygenome.fa -d viruses --data_dir PATH/TO/DATA -o mygenome_virus_result
# Diamond 注释
(eggnog)$ emapper.py -i mygenome.fa -d bact -m diamond --data_dir PATH/TO/DATA -o mygenome_bact_result
```

## 结果

生成的数据文件有：

- \*.emapper.hmm_hits
- \*.emapper.seed_orthologs
- \*.emapper.annotations
