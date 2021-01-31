# Circlize 教程

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

Circlize 是一个生成 circos 类数据可视化的 R 软件包。

```R
# 生成一个随机数，让生成的数据一致
set.seed(12345)

# 生成一个df数据集 data
# letters 是一个常量，包含英文26个小写字母，因此letters[1:8]就是字母a~h
# 取200个这样的字母，replace=T表示可以覆盖，因为8个字母和200个取样数量不同，因此必然又一些字母会重复取样
# sample函数会随机生成符合调价的元素排列
# rnorm(200) 生成200个平均数为0,标准差为1的数字元素
# runif(200) 生成200个0～1之间范围的数
data <- data.frame(factors=sample(letters[1:8], 200, replace=T), x=rnorm(200), y=runif(200))

# 设置全局变量height为0.1，就是每个track占10%宽度
circos.par("track.height" = 0.1)

# 初始化图，类似 ggplot2 的 ggplot 函数
circos.initialize(factors = df$factors, x = df$x)

# 接下来就是一个个track绘制的过程
# 首先生成一个 track，这个 track 包含8个 cell
# circos.track函数其实是trackPlotRegion函数的快捷表达
# y 为 y 轴方向上的数据，x 必须设置了 panel.fun 才能添加
# panel.fun 是用来往cell中添加图示的函数， CELL_META是常量，表示当前cell的meta信息，输入nanmes(CELL_META)可以查看所有CELL_META字段
# uy 是 convert_y 函数的快捷表达
# circos.text用来在涂上标记a~h字符，CELL_META$xcenter标记cell的x位置，第二个元素标记y位置，第三个元素标记字符内容
# circos.axis用来在track上标记ruler标尺
circos.track(factors=data$factors, y=data$y,
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), CELL_META$sector.index)
        circos.axis(labels.cex = 0.6)
	}
)

# 设置一组颜色，给每个cell上色，这里用2个颜色复制4次
col = rep(c("#350E2A", "#9DBC14"), 4)

# 向直接绘制的 track 中的 cell 中加入元素，pch=16为实心圆点，具体参见points函数
# cex表示点的大小，数值越小图示符号越小
circos.trackPoints(data$factors, data$x, data$y, col=col, pch=16, cex=0.5)

# 向之前绘制的 track 中的 cell a 中加入文字
circos.text(-1, 0.5, "text", sector.index="a", track.index=1)

# 设置一组颜色，给新的track上色
bgcol = rep(c("#69EEAD", "#C31705"), 4)

# 数据来源还是data的x值在具体范围内的分布，bin.size设置柱状图的宽度，值越小宽度越小
circos.trackHist(data$factors, data$x, bin.size=0.2, bg.col=bgcol, col=NA)

# 
```

使用 EMBOSS 的 isochore 可以计算sliding windows的GC含量



- circos.points: cell中添加点，比如突变位点，RNA-seq丰度位点
- circos.lines:  cell中添加线，比如GC含量
- circos.segments: adds segments in a cell.
- circos.rect(): adds rectangles in a cell.
- circos.polygon(): adds polygons in a cell.
- circos.text(): adds text in a cell.
- circos.axis() ands circos.yaxis(): add axis in a cell.


- circos.initialize(): allocates sectors on the circle.
- circos.track(): creates plotting regions for cells in one single track.
- circos.update(): updates an existed cell.
- circos.par(): graphic parameters.
- circos.info(): prints general parameters of current circular plot.
- circos.clear(): resets graphic parameters and internal variables.


circos.par 可设置的全局参数

- start.degree 起始位置，默认为0，即时钟3点位置。设置的数值按照时钟逆时针方向旋转，比如一般绘制基因组从12点开始，就设置成 90
