# Perl

---

!!! Abstract "内容简介"
    Perl 是老牌的脚本语言，也是最早的广泛应用在生物信息学软件开发中的语言。虽然近年来像python，R等其他语言越来越多的应用到生物软件中，但是仍然有许多老牌的生信专家在使用perl语言，初步了解基本的perl语法可以帮助理解代码，或者自己开发简单的流程脚本。

## 基本语法

### 1. 数据类型

```Perl
# 变量赋值
$x = 100;
$var = "my data is over";

# 数组赋值
@myarr = (1, 2, 3, 4, 5);
print @myarr[0]; # 1
print @myarr[4]; # 5

# 哈希，相当与python的字典
%dict = ('t1' => 100, 't2' => 300);
print %dict{'t1'}; # 100
```

**字符转义**

| 转义字符 | 含义 |
| -------- | ---- |
| \\ | 反斜线 |
| \' | 单引号 |
| \" | 双引号 |
| \a | 系统响铃 |
| \b | 退格 |
| \f | 换页符 |
| \n | 换行 |
| \r | 回车 |
| \t | 水平制表符 |
| \v | 垂直制表符 |
| \0nn | 创建八进制格式的数字 |
| \xnn | 创建十六进制格式的数字 |
| \cX | 控制字符，x可以是任何字符 |
| \u | 强制下一个字符为大写 |
| \l | 强制下一个字符为小写 |
|\U | 强制将所有字符转换为大写 |
|\L | 强制将所有的字符转换为小写 |
|\Q | 将到\E为止的非单词（non-word）字符加上反斜线 |
|\E | 结束\L、\U、\Q |

### 2. 条件判断

```Perl
# if 结构
if( $i < 100){
  printf "the result is $i";
}

# if else 结构
if( $i < 100){
  printf "the result is $i";
}else{
  printf "wrong number";
}

# if elsif else 结构
if( $i == 100 ){
  printf "a";
}elsif( $i == 200 ){
  printf "b";
}else{
  printf "c";
}

# switch
use Switch;
$var = 10
switch($var){
  case 1 {print "yes"}
  case 5 {print "no"}
  case 10 {print "what?"}
  else {print "wrong"}
}

# 三元表达式
$t = 2000
$current = ($t == 2000)?"Y":"N"; #Y
```

### 3. 循环结构

```Perl
# for 循环
for($i=0;$i<10;$i=$i+1){
  print "loop $i\n";
}

# while 循环
$i = 0
while($i<10){
  printf "i is: $i\n";
  $i = $i + 1;
}

#  do...while 循环
$a = 10;
do{
   printf "a 的值为: $a\n";
   $a = $a + 1;
}while( $a < 15 );

# foreach 循环
@list=(1,2,3,4,5)
foreach $i (@list){
    print "i is: $i\n";
}
```

### 4. 函数使用

```Perl
# 定义函数
sub MyFunc{
  # 调用参数
  $n = @_[0] + @_[1];
  print "hello $n\n";
}
# 调用函数
MyFunc(1,2);
```

### 5. 读写

```Perl
# 打开file.txt
open(DATA, "<file.txt") or die "file.txt 文件无法打开 !";

# 读取文件句柄DATA
while(<DATA>){
   print "$_";
}
# 关闭文件
close(DATA) || die "无法关闭文件";

# 读写标准输入输出
$data=<STDIN>
print "$data";
```

### 6. 正则表达式

```Perl
if($data =~ m/^something/){ ... }

# 匹配终端输入，匹配成功则输出
foreach (<STDIN>){
    chomp;
    if (/perl/){
        print "$_ was matched 'gao'\n";
    }
}
```

## 在线教程

* [Perl Tutorial](http://www.perltutorial.org/)
* [Learn Perl](https://www.learn-perl.org/)
* [Perl 中文教程](https://cn.perlmaven.com/perl-tutorial)
* [Awesome Perl](https://github.com/hachiojipm/awesome-perl)
* [Popular Perl Module](https://github.com/kaxap/arl/blob/master/README-Perl.md)
* [Learn Perl in One Video](https://www.youtube.com/watch?v=WEghIXs8F6c)
* [Mastering Perl Scripting](https://www.youtube.com/watch?v=IoLVCEr207w)
