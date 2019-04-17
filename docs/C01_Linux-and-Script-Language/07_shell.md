# Shell

## 字符串变量处理

字符串变量内置符操作比用外部命令如 expr, awk, sed 等要快，所以能内置符号操作的需求，尽可能用它来完成。

```bash
# 字符串赋值
$ a="R1.fastq.gz"
$ echo $a
R1.fastq.gz

# 给变量赋予默认值
$ echo ${b=100}
100
$ b=20
$ echo ${b=100}
20

$ x=1
$ echo ${y-$x}
1
$ echo ${x-$y}
1
$ echo $y

# 返回字符串长度
$ echo ${#a}
11

# 字符串提取
$ echo ${a:3}    # 从第三个字符串后提取
fastq.gz
$ echo ${a:0:2}  # 从第一个字符串提取2个长度字符
R1
$ echo ${a#R1.}  # 从字符串开头删除R1.
fastq.gz
$ echo ${a%.fastq.gz}  # 从字符串结尾删除匹配.fastq.gz的字符串

# 字符串替换
$ echo ${a/R1/1_1} # 字符串R1替换成1_1, 从头开始替换第一个
1_1.fastq.gz
$ echo ${a//R1/1_1} # 字符串R1替换成1_1, 替换所有
$
```

```bash
$ a1=100
$ a2=10
$ echo ${!a*} # 显示以a开头声明的变量
a1 a2

# 递归循环变量
$ for i in `echo ${!a*}`; do echo $i; done
a1
a2
```

```bash
# 字符串比较
$ s=ATCCGAC
$ if [ $s == 'ATCCGAC' ]; then echo "yes!"; fi
yes!

$ s1=TTCCGAC
$ if [ $t == 'ATCCGAC' ]; then echo "yes!"; else echo "no!"; fi
no!
```

## 无聊的脚本

生成一段特定长度的随机DNA序列

```bash
# 生成1000bp的DNA序列文件
$ for i in {1..1000}; do c=(A T C G); echo ${c[`expr $RANDOM % 4`]} >> 1; done
# 取出文件的\n
$ sed -i ':a;N;$!ba;s/\n//g' 1
# 添加序列名称
$ sed -i '1i\>seq' 1
```
