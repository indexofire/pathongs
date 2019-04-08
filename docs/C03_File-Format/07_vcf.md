# VCF 格式数据文件

VCF: Variant Calling File

体现序列位点变化（点突变、删除、结构变化等）的标准格式文件。文件是压缩并索引过的。最早是为了1000基因组计划而设计的。可以用vcftools从基因组数据中创建。

VCF header 区域

```
##fileformat=VCF4.2

```

描述SNP，INDEL和SV结果的文本文件。在GATK软件中得到最好的支持，当然SAMtools得到的结果也是VCF格式，和GATK的VCF格式有点差别。

## VCF的主体结构
VCF文件分为两部分内容：以“#”开头的注释部分；没有“#”开头的主体部分
3.VCF的10列的意义

1. CHROM ： 参考序列名称
2. POS:variant的位置；如果是INDEL的话，位置是INDEL的第一个碱基位置
3. ID:variant的ID;比如在dbSNP中有该SNP的id，则会在此行给出；若没有，则用’.'表示其为一个novel variant
4. REF:参考序列的碱基
5. ALT:variant的碱基
6. QUAL:Phred格式(Phred_scaled)的质量值，表 示在该位点存在variant的可能性；该值越高，则variant的可能性越大；计算方法：Phred值 = -10 * log (1-p) p为variant存在的概率; 通过计算公式可以看出值为10的表示错误概率为0.1，该位点为variant的概率为90%,qual值与p成正比例
7. FILTER:使用上一个QUAL值来进行过滤的话，是不够的。GATK能使用其它的方法来进行过滤，过滤结果中通过则该值为”PASS”;若variant不可靠，则该项不为”PASS”或”.”。
8. INFO:这一行是variant的详细信息，内容很多，以下再具体详述。
9. FORMAT:variants的格式，例如GT:AD:DP:GQ:PL
10. SAMPLES:各个Sample的值，由BAM文件中的@RG下的SM标签所决定


INFO:

  #DP-read depth：样本在这个位置的reads覆盖度。是一些reads被过滤掉后的覆盖度。DP4:高质量测序碱基，位于REF或者ALT前后
  #QD：通过深度来评估一个变异的可信度。Variant call confidence normalized by depth of sample reads supporting a variant         
  #MQ：表示覆盖序列质量的均方值RMS Mapping Quality
  #FQ：phred值关于所有样本相似的可能性
  #AC，AF 和 AN：AC(Allele Count) 表示该Allele的数目；AF(Allele Frequency) 表示Allele的频率； AN(Allele Number) 表示Allele的总数目。
      对于1个diploid sample而言：则基因型 0/1 表示sample为杂合子，Allele数为1(双倍体的sample在该位点只有1个等位基因发生了突变)，
       Allele的频率为0.5(双倍体的sample在该位点只有50%的等位基因发生了突变)，总的Allele为2； 基因型 1/1 则表示sample为纯合的，Allele数为2，Allele的频率为1，总的Allele为2。
  #MLEAC：Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed
  #MLEAF：Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed
  #BaseQRankSum   比较支持变异的碱基和支持参考基因组的碱基的质量，负值表示支持变异的碱基质量值不及支持参考基因组的，
       正值则相反，支持变异的质量值好于参考基因组的。0表示两者无明显差异。
  #FS  使用F检验来检验测序是否存在链偏好性。链偏好性可能会导致变异等位基因检测出现错误。输出值Phred-scaled p-value，值越大越可能出现链偏好性。
  #InbreedingCoeff    使用似然法检验样本间的近交系数（又或者称为近亲关系）。值越高越可能是近亲繁殖。
  #MQRankSum  比较支持变异的序列和支持参考基因组的序列的质量，负值表示支持变异的碱基质量值不及支持参考基因组的，只针对杂合。
       正值则相反，支持变异的质量值好于参考基因组的。0表示两者无明显差异。实际应用中一般过滤掉较小的负值。
  #BaseCounts   所有样本在变异位点ATCG的数量
  #ClippingRankSum  同前面两个类似，负值表示支持变异的read有更的的hard-clip碱基，正值表示支持参考基因组的的read有更多的hard-clip。0最好，无论是正值还是负值都表示可能可能存在人为偏差。
  #ReadPosRankSum    检测变异位点是否有位置偏好性（是否存在于序列末端，此时往往容易出错）。最佳值为0，表示变异与其在序列上的位置无关。负值表示变异位点更容易在末端出现，正值表示参考基因组中的等位基因更容易在末端出现。
  #ExcessHet   检测这些样本的相关性，与InbreedingCoeff相似，值越大越可能是错误。
  #LikelihoodRankSum  评价支持变异和ref的序列与best hyplotype的匹配性，0为最佳值。负值表示支持变异的read匹配度不及支持ref的匹配度，正值则相反。值越大表示越可能是出现了错误。
  #HaplotypeScore    分数越高越可能出现错误。Higher scores are indicative of regions with bad alignments, typically leading to artifactual SNP and indel calls.
  #SOR：也是一个用来评估是否存在链偏向性的参数，相当于FS的升级版。The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It is an updated form of the Fisher Strand Test that is better at taking into account large amounts of data in high coverage situations. It is used to determine if there is strand bias between forward and reverse strands for the reference or alternate allele. The reported value is ln-scaled.
  #IS：插入缺失或部分插入缺失的reads允许的最大数量
  #G3：ML 评估基因型出现的频率
  #HWE：chi^2基于HWE的测试p值和G3
  #CLR：在受到或者不受限制的情况下基因型出现可能性log值
  #UGT：最可能不受限制的三种基因型结构
  #CGT：最可能受限制三种基因型的结构
  #PV4：四种P值得误差，分别是（strand、baseQ、mapQ、tail distance bias）
  #INDEL：表示该位置的变异是插入缺失
  #PC2：非参考等位基因的phred（变异的可能性）值在两个分组中大小不同
  #PCHI2：后加权chi^2，根据p值来测试两组样本之间的联系
  #QCHI2：Phred scaled PCHI2
  #PR：置换产生的一个较小的PCHI2
  #QBD：Quality by Depth，测序深度对质量的影响
  #RPB：序列的误差位置（Read Position Bias）
  #MDV：样本中高质量非参考序列的最大数目
  #VDB：Variant Distance Bias，RNA序列中过滤人工拼接序列的变异误差范围

FORMAT:
  #AD 和 DP：AD(Allele Depth)为sample中每一种allele的reads覆盖度,在diploid中则是用逗号分割的两个值，
        前者对应ref基因型，后者对应variant基因型； DP（Depth）为sample中该位点的覆盖度。
    #GT：样品的基因型（genotype）。两个数字中间用’/'分 开，这两个数字表示双倍体的sample的基因型。0 表示样品中有ref的allele；
         1表示样品中variant的allele； 2表示有第二个variant的allele。因此： 0/0 表示sample中该位点为纯合的，和ref一致； 0/1 表示sample中该位点为杂合的，有ref和variant两个基因型； 1/1 表示sample中该位点为纯合的，和variant一致。
    #GQ：即第二可能的基因型的PL值，相对于最可能基因型的PL值（其PL=0）而言，大于99时，其信息量已不大，因此大于99的全部赋值99。当GQ值很小时，意味着第二可能基因型与最可能基因型差别不大。
    #GL：三种基因型（RR RA AA）出现的可能性，R表示参考碱基，A表示变异碱基
    #DV：高质量的非参考碱基
    #SP：phred的p值误差线
    #PL：指定的三种基因型的可能性(provieds the likelihoods of the given genotypes)。这三种指定的基因型为(0/0,0/1,1/1)，这三种基因型的概率总和为1。
         和之前不一致，该值越大，表明为该种基因型的可能性越小。 Phred值 = -10 * log (p) p为基因型存在的概率。


4.vcf文件的基因型信息
GT:样品的基因型（genotype）。两个数字中间用’/'分 开，这两个数字表示双倍体的sample的基因型。0 表示样品中有ref的allele； 1 表示样品中variant的allele； 2表示有第二个variant的allele。所以：
0/0表示sample中该位点为纯合位点，和REF的碱基类型一致
0/1表示sample中该位点为杂合突变，有REF和ALT两个基因型（部分碱基和REF碱基类型一致，部分碱基和ALT碱基类型一致）
1/1表示sample中该位点为纯合突变，总体突变类型和ALT碱基类型一致
1/2表示sample中该位点为杂合突变，有ALT1和ALT2两个基因型（部分和ALT1碱基类型一致，部分和ALT2碱基类型一致）
AD和DP:AD(Allele Depth)为sample中每一种allele的reads覆盖度,在diploid（二倍体，或可指代多倍型）中则是用逗号分隔的两个值，前者对应REF基因，后者对应ALT基因型
DP(Depth)为sample中该位点的覆盖度，是所支持的两个AD值（逗号前和逗号后）的加和
例如：
1/1:0,175:175—GT:AD(REF),AD(ALT):DP
0/1:79,96:175
1/2:0,20,56:76
这里的三种类型对应的DP值均是其对应的AD值的加和，1/1的175是0+175，0/1的175是79+96，1/2的76是0+20+56
GQ:基因型的质量值（Genotype Quality）。Phred格式（Phred_scaled）的质量值，表示在该位点该基因型存在的可能性；该值越高，则Genotype的可能性越大；计算方法：Phred值=-10log(1-P)，P为基因型存在的概率。（一般在final.snp.vcf文件中，该值为99，为99时，其可能性最大）
PL:指定的三种基因型的质量值（provieds the likelihoods of the given genotypes）；这三种指定的基因型为（0/0，0/1，1/1），这三种基因型的概率总和为1。该值越大，表明为该种基因型的可能性越小。Phred值=-10log(P)，P为基因型存在的概率。最有可能的genotype的值为0

5. VCF第8列的信息
第8列的信息包括18种，都是以“TAG=Value”，并使用分号分隔的形式，其中很多的注释信息在VCF文件的头部注释中给出，下面对常用的TAG进行解释:
AC，AF和AN
AC(Allele Count) 表示该Allele的数目；AF(Allele Frequency) 表示Allele的频率； AN(Allele Number) 表示Allele的总数目。对于1个diploid sample而言：则基因型 0/1 表示sample为杂合子，Allele数为1(双倍体的sample在该位点只有1个等位基因发生了突变)，Allele的频率为0.5(双倍体的 sample在该位点只有50%的等位基因发生了突变)，总的Allele为2； 基因型 1/1 则表示sample为纯合的，Allele数为2，Allele的频率为1，总的Allele为2
DP（reads覆盖度）
表示reads被过滤后的覆盖度
FS
FisherStrand的缩写，表示使用Fisher’s精确检验来检测strand bias而得到的Fhred格式的p值，该值越小越好；如果该值较大，表示strand bias（正负链偏移）越严重，即所检测到的variants位点上，reads比对到正负义链上的比例不均衡。一般进行filter的时候，推荐保留FS<10~20的variants位点。GATK可设定FS参数。
ReadPosRandSum
Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias.当variants出现在reads尾部的时候，其结果可能不准确。该值用于衡量alternative allele（变异的等位基因）相比于reference allele（参考基因组等位基因），其variant位点是否匹配到reads更靠中部的位置。因此只有基因型是杂合且有一个allele和参考基因组一致的时候，才能计算该值。若该值为正值，表明和alternative allele相当于reference allele，落来reads更靠中部的位置；若该值是负值，则表示alternative allele相比于reference allele落在reads更靠尾部的位置。
进行filter的之后，推荐保留ReadPosRankSum>-1.65~-3.0的variant位点
MQRankSum
该值用于衡量alternative allele上reads的mapping quality与reference allele上reads的mapping quality的差异。若该值是负数值，则表明alternative allele比reference allele的reads mapping quality差。进行filter的时候，推荐保留MQRankSum>-1.65~-3.0的variant位点。
