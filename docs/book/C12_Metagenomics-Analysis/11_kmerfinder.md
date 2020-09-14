# KmerFinder 使用

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    本节简单介绍 [KmerFinder]() 工具的使用。[KmerFinder]() 是DTU大学开发的物种鉴定工具，可以对fastq/fasta数据进行扫描，利用kmer频度的方法与数据比对，鉴定物种。这个工具的在线版本是：https://cge.cbs.dtu.dk/services/KmerFinder/

- [在线版分析工具](https://cge.cbs.dtu.dk/services/KmerFinder/)
- [kmerfinder代码](https://bitbucket.org/genomicepidemiology/kmerfinder)

对于普通用户来说，直接使用[在线版分析工具](https://cge.cbs.dtu.dk/services/KmerFinder/)是最方便的方式。但是如果需要鉴定多个样本或数据集，使用本地分析的方式会更为快速。

## 命令行方式

### 1.安装

```bash
# 建立虚拟环境
$ conda create -n kmerfinder
$ conda activate kmerfinder

# 安装kma工具
(kmerfinder)$ conda install kma

# 下载kmerfinder
(kmerfinder)$ git clone https://bitbucket.org/genomicepidemiology/kmerfinder.git

# 下载数据库
(kmerfinder)$ cd kmerfinder && mkdir database
# 这里设置数据库地址为 /dbs/kmerfinder
(kmerfinder)$ mkdir /dbs/kmerfinder && cd /dbs/kmerfinder
# wget 递归下载20190108版本的所有数据
(kmerfinder)$ wget -r ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/20190108_stable/
# 你也可以访问ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/
# 选择最新的数据作为数据库
# md5验证下载数据
(kmerfinder)$ for i in *.gz; do if [ "$(md5sum $i | awk '{print $1}')" != "$(cat ${i}.md5)" ]; then echo $i "data wrong"; fi; done
# 解压缩数据库并删除压缩包
(kmerfinder)$ tar -zxf *.tar.gz
(kmerfinder)$ rm *.tar.gz
```

```bash
# 利用脚本下载
(kmerfinder)$ git clone https://bitbucket.org/genomicepidemiology/kmerfinder_db.git
(kmerfinder)$ cd kmerfinder_db/
(kmerfinder)$ bash install.sh
```

### 2.运行

使用 python 脚本运行

```bash
# 设置数据库路径
(kmerfinder)$ KmerFinderDB=/dbs/kmerfinder
# 对fastq进行扫描
(kmerfinder)$ python kmerfinder.py -i mydata.fastq.gz -o output -db $KmerFinderDB

# 检索 fasta 数据
(kmerfinder)$ cd KmerFinder
(kmerfinder)$ python3 kmerfinder.py -o result -db /dbs/kmerfinder/bacteria.ATG -i genome.fna
# 生成的结果
(kmerfinder)$ cat result/results.spa

# 检索 fastq 数据
(kmerfinder)$ python3 kmerfinder.py -o result -db /dbs/kmerfinder/bacteria.ATG -i Reads.fastq
# 生成的结果
(kmerfinder)$ cat result/results.spa
```

生成结果文件：

- data.json
- results.spa

## docker 提供在线检索

建立服务器，提供web页面访问，用docker安装配置更为方便。

```bash
$ docker run --rm -it -v $HOME/dbs/kmerfinder -v $(pwd):/workdir kmerfinder -h
```

## 参考资料

1. https://bitbucket.org/genomicepidemiology/kmerfinder


[KmerFinder]: https://bitbucket.org/genomicepidemiology/kmerfinder "KmerFinder"
