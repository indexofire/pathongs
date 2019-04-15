# Linux 软件管理

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

## 3. docker方式

docker 是一种虚拟机，但是其与 virtualbox, vm 等的区别是其启动和调用都是毫秒级的，几乎感觉不到其的加载时间，但又与系统有沙盒模式隔离，保证了安全。docker 最大的好处是直接调用他人写的 DOCKERFILE 就可以下载 image，运行生成 container， DOCKERFILE 解决的安装的过程，除了下载安装速度可能略慢，其他几乎不需要更多的学习曲线。安装完成运行 container 即可使用。

```bash
# Archlinux 安装 docker
$ sudo pacman -S docker

# 使用 DockerFile
$ docker 
```
