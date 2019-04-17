# sed

相比较与awk，sed常用于以行为单位进行替换。

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
