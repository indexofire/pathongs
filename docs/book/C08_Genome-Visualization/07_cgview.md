# CGView 教程

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

cgview是一个生成环状基因组图片的工具。

```bash
# 创建虚拟环境
$ conda create -n cgview
$ conda activate cgview
(cgview)$ conda install cgview
```

```bash
# 添加 Bioperl 路径到环境变量中，默认安装的路径不包含该安装地址
(cgview)$ export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0/

# 生成xml格式配置文件
(cgview)$ cgview_xml_builder.pl -sequence path/to/seq.gbk -output seq.xml
# 生成png格式图片
(cgview)$ cgview -i seq.xml -o seq.png -f png
# 查看图片
(cgview)$ feh seq.png
```

GC content 默认是计算500bp的结果，滑窗移动
