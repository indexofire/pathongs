## bwa 长片段比对

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    本节介绍[Phil Ashton]()等人最近发表的一篇[论文](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3103.html)如何采用 bwa-mem 应用到长片段测序如牛津纳米孔（ONT）。

尽管bwa-mem的第一个版本与PacBio一起工作，但它产生的对齐过于分散而无法使用。我最初认为bwa-mem使用的长精确种子对PacBio读数的~15％错误率不够敏感，但 [Homolog.us](http://homolog.us) 指出BLASR也使用长精确种子。然后我意识到bwa-mem算法也可以使用PacBio数据。随着越来越多有趣的PacBio数据集的出现，我决定尝试一下。

BWA-MEM 算法有两个主要变化，可以更好地支持PacBio数据。首先，我们必须使用放松的评分矩阵，以便Smith-Waterman（SW）可以在有效匹配上给出正分数。在0.7.9和0.7.10中，评分方案是：match = 2，mismatch = -5，gapOpen = -2和gapExt = -1。其次，我添加了一个启发式来过滤初始种子，以减少不成功的种子扩展。对于PacBio读数，bwa-mem在每个种子周围的小窗口中执行SSE2-SW，然后如果SW得分太小（阈值与选项-W成比例）则拒绝种子。这类似于BLAST的X-dropoff启发式算法。除此之外，bwa-mem还实现了间隙修补启发式算法，即使结果对齐不是最优的，它也会尝试将两个共线局部命中与全局对齐连接起来。这种启发式方法有助于让对齐遍历低质量区域，从而减少碎片。通过这些更改，bwa-mem适用于PacBio数据。

由于其较高的错误率，ONT读取带来了新的挑战。单向（1D）读数的初始释放具有高于30％的错误率。 2D读取要好一些，但现在仍然比PacBio有更多错误。 PacBio模式不太合适，需要进一步改进bwa-mem。

ONT特定的更改相对简单。首先，我们使用较短的种子长度和更宽松的阈值-W作为更高错误率的结果。其次，我们根据最近的论文修改得分矩阵匹配= 1，mismatch = -1，gapOpen = -1和gapExt = -1。对于PacBio来说，这个设置也更好。

bwa-mem的ONT模式在很大程度上与 [LAST](http://last.cbrc.jp/) 相当，[LAST](http://last.cbrc.jp/) 是几个团队推荐的映射器。给定相同的评分系统，两个映射器大多数时间生成相同的SW分数。当得分不同时，[LAST](http://last.cbrc.jp/) 往往会成为赢家 - 在这一小部分比赛中得分不同，bwa-mem更有可能错过低质量的命中或无法将部分对齐扩展到正确的位置（我需要走路）通过这些例子来理解为什么会这样）。对于细菌数据，bwa-mem和LAST的速度也差不多。但是，对于人类PacBio读取，bwa-mem的速度要快一些。它更适合人类数据。

[LAST](http://last.cbrc.jp/) 可能是唯一一个有效且准确地工作的映射器，查询序列范围从100bp到100Mbp，无需太多参数调整。这非常令人印象深刻。截至目前，bwa-mem对于超过~10Mbp的查询效果不佳。


```bash
# bwa version > 0.7.10 支持
$ bwa mem -x ont2d ref.fa reads.fq
```

## minimap2
