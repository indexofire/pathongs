# Gubbins 分析

---

## 安装

```bash
# Ubuntu系统自带gubbins
$ sudo apt install gubbins

# 发行版版本冻结在1.x上，如果要使用2.x，比如不想用fastml的可以从源码编译
$ sudo apt install libtool python3-setuptools python3-dev
$ git clone https://github.com/sanger-pathogens/gubbins
$ cd gubbins
$ ./autoreconf -i
$ ./configure
$ make && sudo make install
$ cd python && sudo python3 setup.py install
```

使用 conda 安装到虚拟环境中

```bash
$ conda create -n gubbins
$ conda activate gubbins
(gubbins)$ conda install gubbins=2.4.1
```

## 使用

### 1. 创建 fasta 格式的基因组比对文件


## Reference
