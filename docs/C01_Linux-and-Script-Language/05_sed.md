# sed

sed 全称 Stream Editor，即流编辑器。想对于 awk 处理每一行中各个字段的列特点，sed 常用于以行为单位的编辑操作，比如字符替换等。

!!! note "基本逻辑"
    sed [OPTION]... {script-only-if-no-other-script} [input-file]...

sed 编辑内容可以是标准输入STDIN，也可以是一个文件。

```bash
# 把pattern1替换为pattern2
# 首先搜索pattern1，然后将其替换为pattern2
$ sed s/pattern1/pattern2/g
```

**文件行首插入新行**

```bash
# 在序列test.fas第一行加入序列名称
$ echo "ATCCGAGTTAG" > test.fas
$ sed -i '1i\>seq' test.fas
$ cat test.fa
>seq
ATCCGAGTTAG
```
