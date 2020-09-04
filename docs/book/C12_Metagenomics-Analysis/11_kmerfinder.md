# KmerFinder

---

!!! Abstract "内容简介"
    本节简单介绍 kmerfinder 工具的使用

- [在线版分析工具](https://cge.cbs.dtu.dk/services/KmerFinder/)
- [kmerfinder代码](https://bitbucket.org/genomicepidemiology/kmerfinder)

## 安装

```bash
# 建立虚拟环境
$ conda create -n kmerfinder
$ conda activate kmerfinder

# 安装依赖包
(kmerfinder)$ conda install kma

# 下载kmerfinder
(kmerfinder)$ git clone https://bitbucket.org/genomicepidemiology/kmerfinder.git

# 下载数据库
(kmerfinder)$ cd kmerfinder && mkdir database

```

## 运行

使用 python 脚本运行

```bash
# 检索 fasta 数据
(kmerfinder)$ python3 kmerfinder.py -o . -db $HOME/dbs/kmerfinder/bacteria.ATG -i genome.fna
# 生成的结果
(kmerfinder)$ cat results.spa
```

使用 docker 运行

## 参考资料
