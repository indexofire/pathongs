# Python

---

!!! Abstract "内容简介"
    随着大数据时代的到来，Python 作为一门简单快速的脚本语言，成为了几乎是目前最热门的编程语言。Python 语法见到，书写快速，被越来越多的新进入生物信息学领域的学者作为首选语言进行程序开发。其丰富的第三方库和活跃的社区也是推动 Python 持续发展的动因之一。你可以很方便的搜索到大量的遍布于各个领域的 Python 学习材料。如果要处理生物数据，或者看懂别人程序代码，掌握一点基本的语法知识是非常必要的。本节主要介绍 Python 的基本语法和与生物信息数据处理相关的软件和库。

## 基本语法

### 1. 变量赋值

```python
a = 1
b = "hello"
```

### 2. 数据类型

```python
# 字符串
s = "abc"
s = "中文"

# 数字
x = 1
y = 0.01

# 列表
a = [1, 2, 3]

# 元组
b = (1, 2, 3)

# 字典
c = {'var1': 100, 'var2': 200}
```


### 3. 导入模块

```python
import Bio
from Bio import SeqIO
```

### 4. 循环控制

```python
for i in range(0, 100):
    ...

while True:
    ...
```

### 5. 条件判断

```python
if a == 100 and b is True or c:
    print("wrong")

if a <100:
    print("abc")
elif:
    print("def")
else:
    print("ghi")
```

### 6. 函数定义与调用

```python
def my_func(*args, **kwargs):
    for k, v in kwargs:
        return k, v

# 运行函数
print("hello world")

from math import sqrt
sqrt(2)
```




## 常用模块

对与生物信息或基因组数据常用的学习和分析模块：

- Jupyter Notebook
-

## 参考资料
