# MicroScope

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

![interface](../../assets/C09/10/interface.png)

!!! Abstract "背景简介"
    `LABGeM`生物信息学团队位于`Genoscope`（法国国家测序中心），开发了`MicroScope`，这是一个基于网络的微生物比较基因组分析和手动功能注释平台。

`MicroScope`采用关系数据库模式（PkGDB）存储语法和功能注释管道的结果以及特定基因组的代谢分析，并建立PkGDB与代谢途径数据库（MicroCyc）的关联。`MicroScope`采用Web界面（MaGe）用于帮助评估给定基因组最佳注释所需的所有可用数据。这些数据包括：

- 序列数据（例如，InterPro域预测），
- 关系数据（即同线结果和代谢网络预测），
- 实验数据（给类相关实验），

此外，`MicroScope`提供了许多有用的数据探索功能，例如允许用户执行（复杂）查询，比较基因组研究和代谢分析。它还提供了每个基因组整体统计数据和信息的摘要视图。最后，`MicroScope`可以用作：

- 社区资源，用于对公众可获得的基因组进行比较分析和注释
- 私有资源，因为私有项目的访问权限可以限制为项目负责人定义的一组有限的注释者。



## 参考资料

[MicroScope](https://www.genoscope.cns.fr/agc/microscope)
[说明文档](https://microscope.readthedocs.io/en/latest/index.html)
