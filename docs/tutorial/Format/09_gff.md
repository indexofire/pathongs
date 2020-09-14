# GFF 格式数据

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

GFF（General Feature Format，通用特征格式）基于Sanger GFF2规范。 GFF格式内有九个必填字段，必须以制表符分隔。 如果字段由空格而不是制表符分隔，则轨道将无法正确显示。 有关GFF格式的更多信息，请参阅Sanger的[GFF页面](http://www.sanger.ac.uk/resources/software/gff/)。

GFF格式字段：

1. seqname  - 序列的名称。
2. source  - 生成数据的程序。
3. feature  - 此类功能的名称。
  - CDS
  - start_codon
  - stop_codon
  - exon
4. start  - 序列中要素的起始位置。第一个基数编号为1。
5. end  - 要素的结束位置。序列长度100bp，则结束位置是100。
6. score  - 得分在0到1000之间。如果此注释数据集的轨迹线useScore属性设置为1，则得分值将确定显示此要素的灰度级别（较高的数字=较深的灰色）。如果没有得分值，请输入`.`
7. strand  - `+`，`-`或`.`
8. frame  - 如果是编码外显子/基因，则应为0-2之间的数字，表示第一个碱基的读数帧。如果该特征不是编码外显子，则该值应为“。”。
9. group  - 具有相同组的所有行都链接在一起成为一个项目。

下面是一个 prokka 注释生成的gff部分

```
chr22	TeleGene	enhancer	10000000	10001000	500	+	.	touch1
chr22	TeleGene	promoter	10010000	10010100	900	+	.	touch1
chr22	TeleGene	promoter	10020000	10025000	800	-	.	touch2
```
