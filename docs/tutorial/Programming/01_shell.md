# Shell 脚本编程

{{ git_page_authors }} 更新于：{{ git_revision_date }}

---

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

# 返回字符串长度
$ echo ${\#a}  # 不需要"\"符号，由于教程代码解析问题，加上\保证文档代码解析正确
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
1_1.fastq.gz
```

```bash
# *星号做为通配符，可表示变量名
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

## 循环

### for 循环

循环是常用的功能，往往被用来批量处理数据。循环在终端下直接运行 oneline shell 脚本或者写 shell 程序常会用到。最常用的循环 loop 方式是用 for 来执行。

```bash
# 示例
$ for i in *.fastq.gz; do bwa aln ref.fasta $i > $i.sai; done
```

测序时常用流水号做为样品名称，比如 miseq 测序结果数据常常如下：

```bash
1_S1_R1_L001.fastq.gz   1_S1_R2_L001.fastq.gz
2_S2_R1_L001.fastq.gz   2_S2_R2_L001.fastq.gz
3_S3_R1_L001.fastq.gz   3_S3_R2_L001.fastq.gz
...
12_S12_R1_L001.fastq.gz 12_S12_R2_L001.fastq.gz
```

```bash
# 用 for 循环遍历的写法1-12
$ for i in 1 2 3 4 5 6 7 8 9 10 11 12; do COMMAND; done
$ for i in {1..12}; do COMMAND; done

# 从1-12,间隔3输出，结果i为1,4,7,10
$ for i in {1..12..3}; do COMMAND; done

# 对于同一种文件格式，可以用正则表达来过滤
$ for i in 1.fq 2.fq 3.fq 4.fq; do COMMAND; done
$ for i in *.fq; do COMMAND; done

# 循环常规用法
$ for i in $(COMMAND); do COMMAND; done

# shell 还可以写成类似C语言的循环方式
$ for ((i=1;i<=10;i++)); do COMMAND; done
```

### while 循环

```bash
# while 语法
while command
do
    statement
done
```

具体例子：

```bash
#!/bin/bash
x=1
while [ $x -le 5 ]
do
    echo "Welcome $x times"
    x=$(($x+1))
done
```

### until 循环

```bash
until command
do
    statement
done
```

具体例子：

```bash
declare -i i=10
declare -i sum=0
until ((i>10))
do
    let sum+=i
    let ++i
done
echo $sum
```

## 无聊的脚本练习

**生成一段特定长度的随机DNA序列**

```bash
# 生成1000bp的DNA序列文件
$ for i in {1..1000}; do c=(A T C G); echo ${c[`expr $RANDOM % 4`]} >> 1.fa; done
```

**删除换行符\n**

```bash
# 去除文件的换行符号是常用的命令，比如把一个60字符长度的多行序列变为一行序列。
# sed 流处理的特殊机制，将一行储存在模式空间时会去除\n，返回时再加上\n
# 因此直接用 s/\n//g 匹配是无法去除的
# 这个命令其实是一个多 sed 命令，省略了{}
# :a 定义一个label，用来处理跳转
$ sed -i ':a;N;$!ba;s/\n//g' 1.fa
```

**给fasta序列取名**

```bash
# 利用 sed 添加到第一行，比如给序列命名
$ sed -i '1i\>seq' 1.fa
```
