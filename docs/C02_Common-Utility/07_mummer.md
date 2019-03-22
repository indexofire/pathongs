# Mummer

Mummer 是一个非常快速的细菌基因组比对工具，通过小片段MUMs比对到参考基因组，获得序列的SNP位点差异。因为其快速的特点，有许多细菌分析流程中都整合了Mummer3。Mummer4 尝试整合了一些新功能和多核性能的开放，但是目前还处于测试版，大部分的分析流程也没有更新 Mummer3 到新的版本。

## Mummer3.23

### 安装

```bash
$ conda create -n mummer
$ conda activate mummer
(mummer)$ conda install mummer
(mummer)$ conda install gnuplot=4.6.0
```

## Mummer4

### 安装

```bash
$ conda create -n mummer4
$ conda activate mummer4
(mummer4)$ conda install mummer4
(mummer4)$ conda install gnuplot=4.6.0
```

### 使用 nucmer

```bash
(mummer4)$ nucmer --maxmatch -t 40 -p test ref.fa query.fa
(mummer4)$ cat test.delta
(mummer4)$ mummerplot -f --postscript --large test.delta
```



```bash
(seq)$ seqkit grep -p seq-id Gxw_9143.fas > seq_id.fas
```

## 参考

1. [如何使用MUMmer比对大片段序列](https://vip.biotrainee.com/d/243-%E5%BA%94%E8%AF%A5%E6%98%AF%E6%9C%80%E8%AF%A6%E7%BB%86%E7%9A%84mummer%E4%B8%AD%E6%96%87%E4%BD%BF%E7%94%A8%E8%AF%B4%E6%98%8E)
