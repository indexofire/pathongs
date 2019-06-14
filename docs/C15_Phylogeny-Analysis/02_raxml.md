# RAxML

---

!!! Abstract "内容简介"
    本节介绍最常用的ML进化树构建软件 RAxML 的基本使用方法。


主要参数:

- -m 序列替代模型
  常用碱基模型:  
  * GTRCAT: GTR approximation with optimization of individual per site substitution rates and classiﬁcation of those individual rates into the number of rate categories speciﬁed by -c. This is only a work-around for GTRGAMMA, so be sure not to compare alternative topologies based on their GTRCAT likelihood values.
  * GTRMIX: This option will make RAxML perform a tree inference (search for a good topology) under GTRCAT. When the analysis is ﬁnished RAxML will switch its model to GTRGAMMA and evaluate the ﬁnal tree topology under GTRGAMMA such that it yields stable likelihood values.
  * GTRGAMMA: GTR model of nucleotide substitution with the Γ model of rate heterogeneity. All model parameters are estimated by RAxML. The GTRGAMMA implementation uses 4 discrete rate categories which represents an acceptable trade-off between speed and accuracy. Note that, this has been hard-coded for performance reasons, i.e. the number of discrete rate categories can not be changed by the user.  
  常用氨基酸模型
Available AA models: Values for matrixName: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, GTR. The optional F appendix allows you to  specify whether or not  you want to use empirical base frequencies.
PROTCATmatrixName[F]: (e.g. DAYHOFF or DAYHOFFF) AA matrix speciﬁed by an allowed matrixName, with optimization of individual per–site substitution rates and classiﬁcation of those individual rates into the number of rate categories speciﬁed by -c. This is only a work-around for the GAMMA model of rate heterogeneity, so make sure not to compare alternative topologies based on their PROTCAT-based likelihood values. This is due to the fact that the author assumes that you want to compare trees based on likelihoods if you do a multiple run on the original alignment.
PROTMIXmatrixName[F]: This option will make RAxML perform a tree inference (search for a good topology) under PROTCAT. When the analysis is ﬁnished, RAxML will switch its model to the respective PROTGAMMA model and evaluate the ﬁnal tree topology under PROTGAMMA such that it yields stable likelihood values.
PROTGAMMAmatrixName[F]: AA matrix speciﬁed by matrixName with the  Γ model of rate heterogeneity. All free model parameters are estimated by RAxML. The GAMMA implementation uses 4 discrete rate categories which represents an acceptable trade-off between speed and accuracy. Note that, this has been hard-coded for performance reasons,  i.e. the number of discrete rate categories can not be changed by the user.

-c numberOfCategories

This  option allows you to specify the number of distinct rate categories used into  which the individually optimized rates for each individual site are “thrown”  under -m GTRCAT. The results in [Stamatakis, A.: Phylogenetic models of rate heterogeneity: A high performance computing perspective. In: Proc. of IPDPS2006, Rhodos, Greece (2006)] [pdf] indicate that the default  of -c 25 works ﬁne in most practical cases.


-i initialRearrangementSetting

This allows you to specify an initial rearrangement setting for the initial phase of the search algorithm. If you specify e.g. -i 10 the pruned subtrees will be inserted up to a distance of 10 nodes away from their original pruning point. If you don’t specify -i, a “good” initial rearrangement setting will automatically be determined by RAxML (see Section 5.1 for further details).
We do not support the following sub-options at present:

[-a weightFileName]
[-b bootstrapRandomNumberSeed]
[-d infer a complete random starting tree]
[-e likelihoodEpsilon]
[-f algorithm (d is the default) b|c|d|e|o|s] [-g groupingFileName]
[-j write intermediate trees to a separate file]
[-k optimize branches and model parameters on ootstrapped trees]
[-o outgroupName(s)] [-q multipleModelFileName]
[-r constraintFileName]
[-t userStartingTree]
[-w workingDirectory]
[-v display version information]
[-y ]
[-z multipleTreesFile]
[-# numberOfRuns]

Prior to implementation of the Guigen applet, RAxML was run with only two flags:
-c number of categories; =25; and -m model Of Evolution; =GTRGAMMA
An explanation of the flags is above (the manual [pdf] provides  more details):
