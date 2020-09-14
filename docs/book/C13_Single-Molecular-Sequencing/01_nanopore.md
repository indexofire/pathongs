# Nanopore Minion/Gridon 测序数据分析

更新于: {{ git_revision_date }}

---

!!! Abstract "内容介绍"
    本节以介绍纳米孔测序下机数据的常用分析路线。

Nanopore 测序的优点：

- 实验操作方便
- 数据实时生成，第一时间分析
- 长读长

Nanopore 测序的

Nanopore 测序因为其特点，因此在 basecalling 和 demultiplexing 上有许多软件可以对其操作。不同软件有基于自己的算法，在速度和准确性上有一定的差异，目前许多软件也在不断的发展和改进中。因此与illumina数据不同，在对其数据进行分析前，这2步操作需要有一个基本的认识。

## Basecalling

**High accuracy basecalling**

如果只是用于病原检测，那么可以直接使用 fast basecalling 模式获得碱基，普通高性能电脑基本就可以满足需求。而对于数据质量更高，特别对单个碱基进行分析的工作，最好采用基于 GPU 的高精度 basecalling。

目前 Minion/Gridion 采用的 MinKNOW 软件使用的是 Guppy 进行 basecalling 和 demultiplexing。可以使用 fast 和 accuracy 2种模式，前者在高性能的CPU和SSD上可以达到实时basecalling，后者则需要在一个高性能的GPU和SSD上可以做到实时。

Gridion X5 机器自带了 GV100 GPU 可以实时同时为5张芯片进行basecalling。如果是Minion，则需要外接 GPU 做


## demultiplexing

目前版本的一张 Minion/Gridion 的 Flowcell 上，上样密度优化后可以获得20G以上的碱基数据。虽然芯片可以清洗使用多次，但次数和时间也有限制。实际工作中往往需要在一张芯片同时进行多个样本的测序，这就需要我们在实验室流程中采用连接法或者PCR法把barcode生成到测序序列上，以区分样本。

由于 Nanopore 测序准确度相比 illumina 低很多，因此准确做 demultiplexing 的难度也会大很多。特别是如果使用的是连接酶连接街头，实际获得的连接序列会比较复杂。目前有一些软件可以进行

如果一张芯片上一次进行多个样本的测序，那么就要家barcode序列区分样本，分析时如果拿到的数据没有做过demultiplexing，那么需要自己手工做。这个操作可以有多个软件实现。

官方的 guppy 可以进行 demultiplexing，也可以采用第三方软件对 fast5 进行 demultiplexing 工作。

- deepbinner
- porechop
- qcat
- guppy

deepbinner 与 albacore、porechop 不同，他通过CNN神经网络进行 demultiplexing，目前最新版的MinKnow也使用了Guppy等RNN工具提升basecalling，精度比过去的软件有所提升。

```guppy



```bash
# porechop 使用 -b 参数，可以对reads进行demultiplexing操作
# 解开的各个样本序列放在ouput目录中
(nanopore)$ porechop -i data.fq.gz -b output
```



## 组装分析

### 流程一

**使用软件**

- porechop v0.2.3: 取出接头序列
- NanoFilt v2.2.0: 过滤低质量序列
- canu v1.4: 三代组装软件
- unicycler v0.4.7: 混合组装软件

!!! Warning
    由于 Nanopore 硬件、芯片、软件还在不断改进中，因此第三方分析软件版本和功能可能也会即时变化。

MinKNOW 作为测序软件对 minION 进行测序实验时，测序的默认设置时会做 BASECALLING，生成 fast5 和 fastq 文件。

```bash
# 生成 conda 环境
$ conda create -n nanopore
$ conda activate nanopore
(nanopore)$ conda install porechop nanofilt canu unicycler
```

默认设置时，每4000条通过校验的reads会basecalling成一个fastq文件，这个例子里我们只对12个fastq文件进行分析，每个文件有4K条reads。

```bash
# 去除接头，将fastq文件压缩，应该没有压缩，还是纯文本
# -t 8 为采用8个线程，-i 为输入序列， -o 为生成的去除接头后的fastq文件
(nanopore)$ for fq in {0..11}; do porechop -t 8 -i *_${fq}.fastq -o ${fq}.fastq; done
# 或者将文件合并，在取出接头
# 合并后直接用相应程序对data.fq进行操作，
(nanopore)$ for fq in {0..11}; do cat $fq >> data.fq; done
(nanopore)$ porechop -t 8 -i data.fq -o data_porechop.fq

# 去除低质量碱基，将reads平均Q值低于9的去除
# 随着测序芯片的发展和basecalling的准确度提升，nanopore Q值不断在改进中
# 过滤 -q 参数设置为多少比较合适一般根据你进行的分析来决定
# 对于物种鉴定等reads检测的，一般 -q 7 即可
# 对于精度要求高的溯源，建议 -q 9 以上
(nanopore)$ for trim in {0..11}; do cat ${trim}.fastq | NanoFilt -q 9 >> trimmed.fastq; done

# 采用 canu 直接对纯 nanopore reads 进行组装
# -correct 先进行矫正
# genomeSize 为所拼接物种预估的基因组大小，这里为6M大小
(nanopore)$ canu -correct genomeSize=6m -p vp1 -nanopore-raw trimmed.fastq
# -assemble 组装序列
(nanopore)$ canu -assemble genomeSize=6m -p vp1 -d contigs -nanopore-corrected vp1.correctedReads.fasta.gz

# unicycler 混合组装illumina和nanopore数据
# 这里 illumina PE 测序 reads 为 R1.fastq.tz R2.fastq.gz
(nanopore)$ unicycler -t 8 -1 R1.fastq.gz -2 R2.fastq.gz -l minion.fastq.gz -o hybrid
```

### 目前主流的Nanopore长序列组装软件

- [miniasm](https://github.com/lh3/miniasm): Li heng 开发的长序列组装工具
- [canu](https://github.com/marbl/canu): marbl开发的组装工具
- [wtdbg2](https://github.com/ruanjue/wtdbg2): 国内大牛开发的工具
- [flye](https://github.com/fenderglass/Flye): 最近新出现的热门的工具
- [unicycler](https://github.com/rrwick/Unicycler): 主流的混合组装软件
- [spades](http://cab.spbu.ru/software/spades/): 主流的二代测序组装软件，可支持混合组装

---

## 参考资料

1. https://github.com/rrwick/Porechop
2. https://github.com/rrwick/Deepbinner
