# Linux 基本命令

在 Linux 操作系统命令行下使用生物信息学软件除了软件本身外，还要使用许多 Linux 命令。如果对最基本的 Linux 下如何查看文件，复制粘贴等操作，以 bash 为例。其他的 shell 比如 csh, tcsh 等大部分命令和参数都是一致的。

## 文件操作

- `ls`: 列出文件夹的文件

```bash
# 显示 '~/data' 文件夹里的文件
$ ls ~/data

# 显示包括隐含文件的 '~/data' 文件夹里的所有文件
$ ls -la ~/data

# 文件按照时间排序
$ ls -t ~/data

# 查看所有 fastq 文件
$ ls ~/data/*.fastq
```

- `cd`: 改变文件夹

```bash
# 从当前文件夹 '~/data' 更换到 '~/app' 文件夹
~/data$ cd ~/app

# 从任何位置切换回用户主目录 '/home/User'
$ cd

# 回退上一个访问目录
$ cd -

# 返回上一级目录
$ cd ..

# 返回多级目录冰前往另一个文件夹，这里 '\' 是转义符，将空格键正确转义
$ cd ../../../new\ directory
```

- `pwd`: 显示当前目录

```bash
# 当不知道当前处于什么路径时，可以用这个命令显示
$ pwd
```

- `mkdir`: 建立新文件夹

```bash
# 当前路径下新建一个名叫 'new' 的文件夹
$ mkdir new
```

- `rmdir`: 删除文件夹

```bash
# 删除当前路径的文件夹 'new'
$ rmdir new
```

- `rm`: 删除文件

```bash
# 删除当前文件夹的所有后缀是 '.sra' 的文件
$ rm *.sra

# 删除文件夹 'new' 以及 'new' 下的所有子文件与子目录
$ rm -R new

# 不弹出删除确认提示，删除所有 '.tmp' 文件
$ rm -f *.tmp
```

- `cp`: 复制文件

```bash
# 复制 'test.txt' 文件到文件夹 '~/abc' 中
$ cp test.txt ~/abc
```

- `mv`: 移动文件或文件夹

```bash
# 移动 'test.txt' 文件到文件夹 '~/abc' 中并改名叫 'test1.txt'
$ mv test.txt ~/abc/test1.txt
```

- `which`: 查找可执行文件的系统路径

```bash
# 打印出系统带的 python 程序的路径
$ which python
```

- `wc`: 统计一个文件的行，字符和字节数

```bash
# 输出文件 'text.txt' 的行数，字符数和字节数。
$ wc text.txt
```

- `find`: 查找文件

```bash
# 查找一个文件，比如 .bashrc, /etc路径下的hosts文件
$ find .bashrc
$ find /etc/hosts

# 查找当前目录及子目录后缀为gz的文件
$ find . -type f -name "*.gz"

# 查找子目录深度为4层目录的fasta文件
$ find . -type f -maxdepth 4 -name "*.fasta"

# 查找/var/log路径下权限为755的目录
$ find /var/log -type d -perm 755

# 按照用户查找用户，在tmp路径下查找属于用户 nginx 的文件
$ find /tmp -user nginx
# 用户主home路径下查找属于 root 权限，属性为644的文件
$ find . -user root -perm 644
```

- `ps`: 查看系统进程

```bash
#ps会在终端打印系统进程，各列的含义是:

* USER: 运行该进程的用户
* PID: 运行着的命令(CMD)的进程编号
* %CPU: CPU占用
* %MEM: 内存占用
* VSC:
* RSS:
* TTY: 命令所运行的位置（终端）
* STAT:
* TIME: 运行着的该命令所占用的CPU处理时间
* COMMAND: 该进程所运行的命令
```

```bash
# 显示详细的进程信息
$ ps -waux

# 过滤用户root的进程
$ ps -u root

# 根据不同参数使用来排序进程，并只现实排名前10的进程
$ ps -aux --sort -pcpu | head -n 11
$ ps -aux --sort -pmem | head -n 11
$ ps -aux --sort -pcpu,+pmem | head -n 11

# 过滤进程名
$ ps -f -C chrome

# 根据PID过滤
$ ps -L 1000

# 树形现实进程
$ ps -axjf
```

## 输出文件头尾部内容

操作大文件如fastq测序原始数据，常常不需要将几百M甚至上G的数据全部写入到内存打开文件，只需要了解头部数据即可。用head就非常轻便。head默认输出前10行内容。

```bash
# 输出文件前5行内容
$ head -5 text.txt
```

- `tail`: 输出文件尾部内容

```bash
# 输出文件的最后5行内容
$ tail -5 text.txt
```

- `cat`: 输出文件内容

```bash
# 显示文件 'text.txt' 内容
$ cat text.txt
```

- `nl`: 将文件标行

```bash
$ cat text.txt | nl
```

- `grep`: 截取输入字符的选定 pattern 并输出所在的行

虽然很多发行版用egrep代替grep，但是其对正则表达式的支持还是不够优秀，所以常用awk代替grep来进行模式匹配。

```bash
# 显示文件 'text.txt' 中含有字符 'abc' 的行
$ cat text.txt | grep 'abc'
$ awk '/abc/' text.txt
```

## 压缩/解压缩

```bash
# 压缩文件 foo 或者文件夹 foo 成为一个 foo.tar.gz 压缩包
$ tar -zcf foo.tar.gz foo

# 解压缩文件 foo.tar.gz
$ tar -zxf foo.tar.gz
```

## 客户端与服务器端交互

1. ssh 登录服务器

```bash
# server_ip为所要登录服务器ip
$ ssh user@server_ip
$ ssh mark@10.44.35.122
```

2. 客户端与服务器端拷贝文件

这里以客户端操作为例，服务器ip假设为10.44.35.122

```bash
# 从服务器拷贝文件 /data/1.fasta 到本地
$ scp mark@10.44.35.122:/data/1.fasta .
# 从本地上传文件 2.fasta 到服务器/data
$ scp 2.fasta mark@10.44.35.122:/data/2.fasta
```

简单的传输单个文件用scp还算方便，但如果频繁交互，或者需要可视化界面来操作，我们需要更方便的工具，这里推荐 sshfs

```bash
# 挂载服务器 /data 到本地linux机器的 ~/server 目录。
$ sshfs mark@10.44.35.122:/data ~/server -o idmap=user -o allow_other
```

对于win平台用户，简单的ssh登录工具常用的如putty，如果要上下传文件，也可以搭建ftp服务来实现。

3. rsync

rsync命令是一个远程数据同步工具，可快速同步多台主机间的文件。rsync使用所谓的“rsync算法”来使本地和远程两个主机之间的文件达到同步，这个算法只传送两个文件的不同部分，而不是每次都整份传送，因此速度相当快。

```bash
$ rsync -avz user@server_ip:/data/* /data
$ rsync -av ~/test.py user@server_ip::/home/user
```

## Reference

1. http://linoxide.com/how-tos/linux-ps-command-examples/
2. https://commandlinefu.cn/
