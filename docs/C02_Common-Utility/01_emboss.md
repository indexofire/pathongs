# EMBOSS

EMBOSS 是欧洲分子生物学组织开发的 Unix/Linux 下的生物学分析工具。EMBOSS 包含工具众多，这里只介绍与微生物基因组分析可能会用到的一些工具，所有的软件和其文档参考官方文档。

[软件列表](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/)

## 安装 EMBOSS

```bash
$ conda create -n emboss emboss
$ conda activate emboss
```

EMBOSS 工具集包含众多基于命令行的工具，可以集成到分析工具流中。

### 序列处理

**seqret**: 生成序列或格式化序列

```bash
# 读入序列，将其格式化为规范的fasta格式数据
(emboss)$ seqret -sequence 1.fas -outseq 2.fas

# 通过交互方式进行操作
(emboss)$ seqret
Read and write (return) sequences
Input (gapped) sequence(s):         <--- 输入序列文件名称
output sequence(s) [...]:           <--- 输入序列保存的文件名称
```

**transeq**: DNA/RNA序列翻译成氨基酸序列，文件以 `.pep` 后缀保存  
**backtranseq**: 氨基酸序列转换成DNA序列。

```bash
# 将 1.fas 序列翻译成氨基酸序列保存成 1.pep
(emboss)$ transeq 1.fas 1.pep
# 将 1.pep 氨基酸序列转换成碱基格式保存成 1.fas
(emboss)$ backtranseq 1.pep 1.fas
```

### 引物设计

**eprimer3**:
**primersearch**

### 双序列比对

双序列(Pairwise Alignment)比对的软件

**needle**:
**water**:

```bash
(emboss)$ water seq1.fas seq2.fas seq1v2.alignment
```


**needleall**

### 多序列比对

**emma**:
**edialign**:
