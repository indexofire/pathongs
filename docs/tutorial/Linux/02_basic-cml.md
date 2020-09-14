# Linux 基本命令

{{ git_page_authors }} 更新于: {{ git_revision_date }}

![Linux基本命令](../assets/images/C01/01/banner.jpg)

---

!!! Abstract "简介"
    在 Linux 操作系统命令行下使用生物信息学软件除了需专门的分析软件外，往往还要结合使用许多 Linux 命令；而且许多 Linux 命令可以通过管道符号的方式构建数据处理流程。因此对最基本的 Linux 文件查看、复制、粘贴等基本操作必须应掌握。

对于没有只接触过 windows 图形界面的用户，需要知道几个基本知识点：

!!! info "知识点"
    1. 路径: 和 windows 盘符概念不同，Linux 只有类似树的根状结构一般的目录管理，所有路径都挂在`/`下。和 windows 路径分割符号`\`不同，Linux 是使用`/`。路径中'.'表示当前目录，'..'表示上一级目录。路径分为绝对路径和相关路径。
    2. 用户: Linux 系统多用户管理的权限区别，windows 下个人用户一般都使用管理员权限直接运行，而 Linux 下应严禁使用 root 用户日常操作。所以一般都给一个 wheel 用户组用户以 sudo 权限，以便必要是提供 sudo 来进行管理员操作。
    3. 后缀: windows 下习惯用3位字符来区别不同文件类型。Linux 没有这个必要，后缀可以任何方式表示，也可以不代表任何意义。文件一般只有文本文件和二进制文件2种。前者可以用文本编辑器打开编辑。
    4. 发行版: Linux 有数量众多的发行版，虽然内核使用的都是 Linux，但在软件管理等方式方面往往使用的是不同方式，所以会根据发行版有不同命令来安装软件或更新系统。但对于本节介绍的基本命令，是与发行版无关的，任何桌面 Linux 系统都具有的命令。

## 文件操作

- `ls`: 列出目录
- `cd`: 改变目录
- `pwd`: 显示当前目录
- `mkdir`: 建立新目录
- `rmdir`: 删除目录
- `rm`: 删除文件
- `cp`: 复制文件
- `mv`: 移动文件或目录
- `which`: 查找可执行文件的系统路径
- `wc`: 统计一个文件的行，字符和字节数

```bash tab="ls"
# 显示 '~/data' 目录里的文件
$ ls ~/data

# 列表方式显示包括隐含文件的 '~/data' 目录里的所有文件
$ ls -la ~/data

# 文件按照时间排序
$ ls -t ~/data

# 查看所有后缀为 fastq 的文件
$ ls ~/data/*.fastq
```

```bash tab="cd"
# 从当前目录切换到 '~/app' 目录
$ cd ~/app

# 从任何位置切换回用户主目录 '/home/User'
$ cd

# 回退上一个访问目录
$ cd -

# 返回上一级目录
$ cd ..

# 使用相对路径，进入另一个目录"new directory"，这里 '\' 是转义符，将空格键正确转义
$ cd ../../../new\ directory
```

```bash tab="pwd"
# 当不知道当前处于什么路径时，可以用这个命令显示
$ pwd

# 需要使用当前路径到命令中，如建立links
$ ln -s `pwd`/mylib.so /usr/lib/
```

```bash tab="mkdir"
# 当前路径下新建一个名叫 'new' 的目录
$ mkdir new
```

```bash tab="rmdir"
# 删除当前路径的目录 'new'
$ rmdir new
```

```bash tab="rm"
# 删除当前目录的所有后缀是 '.sra' 的文件
$ rm *.sra

# -R 参数删除目录
# 删除目录 'new' 以及 'new' 下的所有子文件与子目录
$ rm -R new

# 不弹出删除确认提示，删除所有 '.tmp' 文件
$ rm -f *.tmp

# 千万不能做的事情！也是linux频道里经常能看到的笑话
$ sudo rm -rf /*
```

```bash tab="cp"
# 复制 'test.txt' 文件到目录 '~/abc' 中
$ cp test.txt ~/abc

# -R 复制目录
# 复制 mydir 到用户目录下
$ cp -R mydir/ ~
```

```bash tab="mv"
# 移动 'test.txt' 文件到目录 '~/abc' 中并改名叫 'test1.txt'
$ mv test.txt ~/abc/test1.txt
```

```bash tab="which"
# 打印出系统带的 python 程序的路径
$ which python
```

## 系统查看

- `find`: 查找文件
- `ps`: 查看系统进程
- `top`: 查看当前系统状态和进程


```bash tab="find"
# 查找一个文件，比如 .bashrc, /etc路径下的hosts文件
$ find .bashrc
$ find /etc/hosts

# 查找当前目录及子目录后缀为gz的文件
$ find . -type f -name "*.gz"
# 查找所有的 html 网页文件
$ find . -type f -regex ".*html$"
# 将所有 fasta 文件中序列名字包含node的改成seq
$ find . -type f -regex ".*fasta$" | xargs sed 's/NODE/Seq/g'

# 查找子目录深度为4层目录的fasta文件
$ find . -type f -maxdepth 4 -name "*.fasta"

# 查找/var/log路径下权限为755的目录
$ find /var/log -type d -perm 755

# 按照用户查找用户，在tmp路径下查找属于用户 nginx 的文件
$ find /tmp -user nginx
# 用户主home路径下查找属于 root 权限，属性为644的文件
$ find . -user root -perm 644
```

```bash tab="ps"
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

```bash tab="top"
$ top
```

!!! info "进程列表"

    ps会在终端打印系统进程，各列的含义是:

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



## 显示文本文件内容

操作大文件如fastq测序原始数据，常常不需要将几百M甚至上G的数据全部写入到内存打开文件，只需要了解头部数据即可。用head就非常轻便。head默认输出前10行内容。

- `head`
- `tail`
- `cat`
- `less`
- `more`
- `wc`:
- `nl`:
- `grep`: 截取输入字符的选定 pattern 并输出所在的行

```bash tab="head"
# 输出文件前5行内容
$ head -5 text.txt
```

```bash tab="tail"
# 输出文件的最后5行内容
$ tail -5 text.txt
```

```bash tab="cat"
# 显示文件 'text.txt' 内容
$ cat text.txt
```

```bash tab="less"
# 显示文件 'text.txt' 内容，支持分页
$ less text.txt
```

```bash tab="more"
$ more text.txt
```

```bash tab="wc"
# 输出文件 'text.txt' 的行数，字符数和字节数。
$ wc text.txt

# 输出行数
$ wc -l text.txt
```

```bash tab="nl"
# 输出行号
$ cat text.txt | nl
```

```bash tab="grep"
# 显示文件 'text.txt' 中含有字符 'abc' 的行
$ cat text.txt | grep 'abc'

# 输出文件中行以>开头的行
$ grep '^>' myseq.fasta
```

## 压缩/解压缩文件

对序列数据常需要进行压缩处理，可以使用 gzip/tar， 特殊格式压缩可以使用第三方软件，个人常用7z来进行。

```bash
# gzip 压缩文件 1 为 1.gz
$ gzip 1
# gunzip 解压缩
$ gunzip 1.gz

# 压缩文件 foo 或者目录 foo 成为一个 foo.tar.gz 压缩包
$ tar -zcf foo.tar.gz foo
# 解压缩文件 foo.tar.gz
$ tar -zxf foo.tar.gz
```

## 客户端与服务器端交互

Linux 服务器上计算获得的结果，可以通过以下常用方式方便的下载和上传。

1. ssh 登录服务器

安装了 openssh 后，直接通过 ssh 登录开启了 sshd 服务的 Linux 服务器，直接在远端进行操作。

```bash
# server_ip为所要登录服务器ip
$ ssh user@server_ip
```

2. 客户端与服务器端拷贝文件

这里以客户端操作为例，服务器ip假设为10.44.35.122

```bash
# 从服务器拷贝文件 /data/1.fasta 到本地
$ scp mark@10.44.35.122:/data/1.fasta .
# 从本地上传文件 2.fasta 到服务器/data
$ scp 2.fasta mark@10.44.35.122:/data/2.fasta
```

rsync命令是一个远程数据同步工具，可快速同步多台主机间的文件。rsync使用所谓的“rsync算法”来使本地和远程两个主机之间的文件达到同步，这个算法只传送两个文件的不同部分，而不是每次都整份传送，因此速度相当快。

```bash
# 传送 server_ip 服务器上的文件到本地 Linux 系统 /data 路径中
$ rsync -avz user@server_ip:/data/* /data
# 传送本地 test.py 文件到服务器 server_ip 用户目录中
$ rsync -av ~/test.py user@server_ip::/home/user
```

3. 简单的传输单个文件用scp还算方便，但如果频繁交互，或者需要可视化界面来操作，我们需要更方便的工具，对于win平台用户，简单的ssh登录工具常用的如putty，如果要上下传文件，也可以搭建ftp服务来实现。这里推荐 sshfs。

```bash
# 挂载服务器 /data 到本地linux机器的 ~/server 目录。
$ sshfs mark@10.44.35.122:/data ~/server -o idmap=user -o allow_other
```

## Reference

1. http://linoxide.com/how-tos/linux-ps-command-examples/
2. https://commandlinefu.cn/
