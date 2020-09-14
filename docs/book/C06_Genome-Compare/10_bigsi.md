# BIGSI

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"

    而 2019年初发表在 Nature Biotechnology 上的`BIGSI`通过建立独特的索引 index，可以实现超高速的 SRA 搜索，很可能未来成为一个很好的搜索细菌 DNA 序列解决方案。

## 1. SRA 数据搜索

有时候从事细菌学研究的我们都想要一个搜索 DNA 的 google 工具，只要把想要搜索的DNA序列，比如耐药基因序列复制到剪切框中点击确认就能帮我们找到包含这些耐药基因的数据。过去常常使用 NCBI Blast 来实现，但是它的 nt/nr 数据库无法支持 SRA 数据。blast 算法对于海量的短序列比对效率很低，只能先组装后再纳入数据库。

随着高通量测序数据的大规模提交，我们希望能从这些基因组 SRA 测序数据中去挖掘包含的序列。现有的方法是采用类似 bwt mapping 的比对的方式，这给计算带来了非常高的的要求。

## 2. 安装bigsi

如果从源代码编译安装，会遇到序列依赖问题，推荐使用 bioconda 包直接安装。

```bash
# 创建虚拟环境
$ conda create -n bigsi
$ conda activate bigsi

# 安装 cython
# 另外建议使用 python3.6
(bigsi)$ conda install cython

# 安装 berkeley-db
(bigsi)$ cd $CONDA_PREFIX/share
(bigsi)$ wget http://download.oracle.com/berkeley-db/db-4.8.30.tar.gz
(bigsi)$ tar -zxf db-4.8.30.tar.gz
(bigsi)$ cd db-4.8.30/build_unix
(bigsi)$ ../dist/configure --prefix $CONDA_PREFIX
(bigsi)$ make && make install
(bigsi)$ rm ../../db-4.8.30-tar.gz
# 安装 berkeley-db python 支持依赖库 bsddb3
(bigsi)$ conda install bsddb3

# 安装 rocksdb
(bigsi)$ conda install -c activisiongamescience rocksdb python-rocksdb

# 安装 bigsi
(bigsi)$ conda install bigsi

# 安装 mccortex
(bigsi)$ conda install mccortex

# 安装 mykrobe
# bigsi >= 3.5 可以用 mykrobe 做 variant search
(bigsi)$ conda install mykrobe
# mykrobe 安装脚本下载数据库，需要翻墙，因此可以用 proxychains 工具代理下载
```

## 3. 使用bigsi

### 3.1 从自己的数据开始

首先使用 kmer 工具构建索引，这里使用 mccortex。然后配置数据库，这里使用 berkeleydb。bigsi目前还支持rocksdb 和 redis 等。

```bash
# 建立 kmer
(bigsi)$ mccortex31 build -k 31 -s test1 -1 1_R1.fq.gz test1.ctx
(bigsi)$ mccortex31 build -k 31 -s test2 -1 1_R2.fq.gz test2.ctx
# 选择 berkeleydb 做为 db_backend
(bigsi)$ export BIGSI_CONFIG=berkeleydb.yaml
# 生成 bloom
(bigsi)$ bigsi bloom test1.ctx test1.bloom
(bigsi)$ bigsi bloom test2.ctx test2.bloom
# 合并 bloom
(bigsi)$ bigsi build test1.bloom test2.bloom -s s1 -s s2
# seq 长度要超过 kmer 长度
(bigsi)$ bigsi search CGGCGAGGAAGCGTTAAATCTCTTTCTGACGGGACC
```

### 3.2 构建公共数据库

如果不仅仅想在自己的测序数据集中搜索，而是想到公共数据库中去查询，可以用目前开放的 bigsi server 服务。如果想自己建立服务器，那么可以直接下载已经建立好的细菌索引文件。

```bash
# 下载官方建立的2018年细菌基因组数据
(bigsi)$ wget -c -e robots=off --cut-dirs 6 -m -np http://ftp.ebi.ac.uk/pub/software/bigsi/nat_biotech_2018/all-microbial-index-v03/
# 合并数据
(bigsi)$ cat all-microbial-bigsi-v03* > combined-index

# berkeleydb 设置 berkeleydb.yaml
h: 3
m: 25000000
nproc: 4
k: 31
storage-engine: berkeleydb
storage-config:
  filename: combined-index
  flag: "c"
```

## 参考资料

- https://github.com/Phelimb/BIGSI
