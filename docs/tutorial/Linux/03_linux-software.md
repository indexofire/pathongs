# Linux 软件管理

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

在个人Linux上运行和管理相对比较简单。但在服务器上要想做好多用户管理则要复杂许多。如果要运行一些负载均衡或者运行如Galaxy实例，则需要管理员具备相当丰富的经验和技术。这些技术大部分都是与计算机相关，本节只简单讨论一下与生物软件运行相关的Linux系统管理。

## 1. 系统级安装

Linux 使用比 Windows 非常便捷的一个地方，在于其软件管理（虽然 Windows 目前也有第三方软件如`360软件管家`可以进行自动安装与更新）。大部分主流发行版都内带了命令行软件管理工具。比如ubuntu的apt工具，archlinux的pacman工具等。Linux 安装系统自带的预编译包的好处是软件仓库自动解决了软件/库依赖的问题。你安装软件A，会自动安装所需要的其他库。

```bash
# Ubuntu 类发行版软件管理
# 软件更新
$ sudo apt update
$ sudo apt upgrade
# 软件安装
$ sudo apt install pkg-name
# 软件删除
$ sudo apt remove pkg-name
# 软件搜索
$ sudo apt search pkg-name

# Archlinux 类发行版软件管理
# 软件更新
$ sudo pacman -Syu
# 软件安装
$ sudo pacman -S pkg-name
# 软件删除
$ sudo pacman -R pkg-name
# 软件搜索
$ sudo pacman -Ss pkg-name
```

安装系统级或者用户级软件

## 2. conda 虚拟环境

由于很多开发者需要在不同环境下调用不同版本，或者使用不同版本的 java, python进行开发和运行程序。而各种版本的依赖库由于对版本的支持度，会造成错综复杂的结果。因此我们常常将一个应用放在与系统级程序隔离的虚拟环境中。目前最流行的就是使用 conda 工具。 conda 包含 anaconda 和 miniconda 2个版本。前者自带了很多数据分析相关的应用程序和包，后者则是最小化了第三方程序，只提供了一个mini的conda支持。

```bash
# 下载安装 miniconda
$ cd /tmp
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# 只能在当前用户安装 conda，不能安装到系统级权限下
$ sh Miniconda3-latest-Linux-x86

# 设置路径，个人一般设置为 $HOME/.conda
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

建立 conda 虚拟环境

```bash
# 创建虚拟环境 myenv
$ conda create -n myenv
$ conda activate myenv
# 在 myenv 环境中安装软件
(myenv)$ conda install mysoftware
```

## 3. Git 代码管理工具

git 是 Linux 内核之父 Linus 开发的代码管理工具。目前社区中比较活跃的 git 在线社群服务有 github, gitlab, bitbucket(过去是hg,现在也转为git)，国内的例如coding.net,oschina。github是最早最成功的 git 在线管理社区，极好的粘合了广大“码农”，众多主流的开源软件都驻扎在 github 中。

许多生物学软件的源代码也选择托管在 github 中（老牌的可能在 sourceforge.net），如果我们要研究代码，或者使用脚本程序，或者通过源代码安装编译安装软件等，就要学会最基本的 git 命令。而且 github 上还有许多基于 markdown, rst, jupyter, pdf 等格式的教学资料，是学习的最好场所之一。

```bash
# 安装 git
# ubuntu
$ sudo apt install git-core
# archlinux
$ sudo pacman -S git

# 克隆 bwa 源码仓库到本地
$ git clone https://github.com/lh3/bwa.git
$ cd bwa

# 编译安装 bwa
$ make
$ sudo cp bwa /usr/local/sbin

# 使用 bwa
$ bwa mem ...
```

## 4. docker方式

docker 是一种虚拟机，但是其与 virtualbox, vm 等的区别是其启动和调用都是毫秒级的，几乎感觉不到其的加载时间，但又与系统有沙盒模式隔离，保证了安全。docker 最大的好处是直接调用他人写的 DOCKERFILE 就可以下载 image，运行生成 container， DOCKERFILE 解决的安装的过程，除了下载安装速度可能略慢，其他几乎不需要更多的学习曲线。安装完成运行 container 即可使用。

**Ubuntu 18.04 Server安装docker**

```bash
# 添加 docker 仓库
$ sudo apt update
$ sudo apt install apt-transport-https ca-certificates curl software-properties-common
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
$ sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable"
$ sudo apt update
$ apt-cache policy docker-ce
# 安装 docker
$ sudo apt install docker-ce
# 查看 docker 启动状态
$ sudo systemctl status docker
# 如果没有启动系统是运行 docker

# 将当前用户加入 docker 组，运行 docker 时可以不用输入 sudo 命令
$ sudo usermod -aG docker ${USER}
# 更新组信息
$ su - ${USER}
# 确认当前用户已经加入 docker 组
$ id -nG
```

**docker镜像下载**

```bash
# 运行软件自带的 hello-world 镜像
$ docker run hello-world
# 镜像搜索
$ docker search ubuntu
# 镜像下载
$ docker pull ubuntu
# 查看当前已有的镜像
$ docker images
```

**docker容器运行**

```bash
# 以交互方式运行容器 ubuntu
$ docker run -it ubuntu

# 查看当前运行容器
$ docker ps
# 查看所有容器
$ docker ps -a
```
