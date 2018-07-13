# Mummer

## 安装 Mummer4

```bash
$ conda create -n mummer4
$ conda activate mummer4
(mummer4)$ conda install mummer4
(mummer4)$ conda install gnuplot=4.6.0
```

## 使用 Mummer4

### 使用 nucmer

```bash
(mummer4)$ nucmer --maxmatch -t 40 -p test ref.fa query.fa
(mummer4)$ cat test.delta
(mummer4)$ mummerplot -f --postscript --large test.delta
```

```bash
(seq)$ seqkit grep -p seq-id Gxw_9143.fas > seq_id.fas
```


## 安装 Mummer3

```bash
$ conda create -n mummer
$ conda activate mummer
(mummer)$ conda install mummer
(mummer)$ conda install gnuplot=4.6.0
```



## 参考

1. [如何使用MUMmer比对大片段序列](https://vip.biotrainee.com/d/243-%E5%BA%94%E8%AF%A5%E6%98%AF%E6%9C%80%E8%AF%A6%E7%BB%86%E7%9A%84mummer%E4%B8%AD%E6%96%87%E4%BD%BF%E7%94%A8%E8%AF%B4%E6%98%8E)
