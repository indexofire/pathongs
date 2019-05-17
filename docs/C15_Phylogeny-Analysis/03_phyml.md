# PhyML 使用

[PhyML](https://github.com/stephaneguindon/phyml/) 是由 Stephane Guindon 开发的一款进化树软件。

## 1. 安装

### 1.1 从源代码安装

```bash
$ git clone https://github.com/stephaneguindon/phyml
$ cd phyml && sh ./autogen.sh
# 安装 phyml 以及相关的软件 phyrex 和 phytime
$ ./configure --enable-phyml
$ make
$ sudo make install
```

### 1.2 预编译包安装

```bash
$ conda activate phylo
(phylo)$ conda install phyml
```

## 2. 使用

```bash
# 对碱基序列构建GTR进化树
$ phyml -i nucl.maf -d nt -b 1000 -m GTR -f m -v e -a e -o tlr
# 对氨基酸序列构建JTT进化树
$ phyml -i prot.maf -d aa -m JTT -a e
```

!!! note "软件参数"
    - -i seq_file_name 输入文件，phylip 格式的多序列比对结果。
    - -d data_type default：nt 该参数的值为 nt, aa 或 generic。
    - -b int 设置 bootstrap 次数。
    - -m model 设置替代模型。 核酸的模型有： HKY85（默认的）, JC69, K80, F81, TN93, GTR ; 氨基酸的模型有：LG （默认的）, WAG, JTT, MtREV, Dayhoff, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, HIVw, HIVb 。
    - -f e,m or fA,fC,fG,fT 设置频率计算的方法。 e 表示使用比对结果中不同氨基酸或碱基出现的频率来计算； m 表示使用最大似然法计算碱基频率，或使用替换模型计算氨基酸频率； fA,fC,fG,fT 则是 4 个浮点数，表示 4 中碱基的频率，仅适合核酸序列。
    - -v prop_invar 设置不变位点的比例，是一个[0,1]区间的值。或者使用 e 表示程序获得其最大似然估计值。
    - -a gamma gamma 分布的参数。此参数值是个正数，或者使用 e 表示程序获得其最大似然估计值。在 ProtTest 软件给出的最优模型中含有 G 时，使用该参数。
    - -o params 参数优化的选项。t 表示对 tree topology 进行优化； l 表示对 branch length 进行优化； r 表示对 rate parameters 优化。params=tlr 这表示对 3 者都进行优化。 params=n 表示不进行优化。

## 参考资料

[PhyML Paper](http://www.atgc-montpellier.fr/download/papers/phyml_2010.pdf)
