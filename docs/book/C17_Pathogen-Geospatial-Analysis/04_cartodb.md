# CartoDB

---

[CartoDB](http://carto.com) 是一个地理数据可视化软件。carto.com 提供了基于 Web 的实例化应用，其免费帐号数据使用有一些权限限制，同时制作的 map 必须公开。如果要做私有化数据的 map，或者自己团队内部数据分析分享，就要自己建立服务。而CartoDB是开源的，很方便的建立相关应用和服务。

## 1. 安装

cartodb 作为一套应用包，如果自己手动一个个安装，对于普通用户来说学习曲线是比较高的。因此我们一般使用 docker 来安装和部署服务。

### 1.1 Ubuntu 安装 docker

```bash
# 删除旧版本 docker
$ sudo apt remove docker docker-engine
$ sudo apt purge docker docker-engine

# 安装依赖，如果已经安装可以跳过
$ sudo apt install -y apt-transport-https ca-certificates curl software-properties-common

# 添加官方 GPG key
$ curl -fsSL https://mirrors.aliyun.com/docker-ce/linux/ubuntu/gpg | sudo apt-key add -
$ sudo apt-key fingerprint 0EBFCD88
$ sudo add-apt-repository "deb [arch=amd64] https://mirrors.aliyun.com/docker-ce/linux/ubuntu \
> $(lsb_release -cs) stable"

# 安装 docker
$ sudo apt update
$ sudo apt install -y docker-ce

# 测试安装是否成功，如果看到下载远程镜像信息，并由打印了欢迎信息，则表示安装成功。
$ docker run hello-world

# 用 daocloud 加速镜像下载，并重启服务
$ curl -sSL https://get.daocloud.io/daotools/set_mirror.sh | sh -s http://8ad7943c.m.daocloud.io
$ sudo systemctl restart docker.service

# 添加用户到 docker 组
$ grep 'docker' /etc/group | grep $USER
# 如果没有输出，则可以把当前用户添加到 docker 组，以免运行 docker 时需要 root 权限。
$ sudo usermod -aG docker $USER
```

### 1.2 安装 cartodb

```bash
# 克隆 docker-cartodb 到本地仓库，构建镜像
$ git clone https://github.com/sverhoeven/docker-cartodb.git
$ docker build -t=sverhoeven/cartodb docker-cartodb/
# 接下来就是漫长的等待完成虚拟系统的安装过程

# 运行镜像，ip_address输入服务器ip地址，让docker container绑定到服务器80端口上，用浏览器访问服务器ip即可看到carto web界面
$ docker run -d -p 80:80 -h ip_address sverhoeven/cartodb
```

## 2. 使用

```bash
# docker数据层是独立的，所以可以建立单独的PostgreSQL persistant
$ docker create --name cartodb_pgdata sverhoeven/cartodb

# Change to directory to save the Postgresql data dir (cartodb_pgdata) of the CartoDB image
$ docker cp cartodb_pgdata:/var/lib/postgresql $PWD/cartodb_pgdata
$ docker rm -f cartodb_pgdata

# Inside container cartodb_pgdata is owned by postgres (uid=105) user,
# it should be owned by same user on the local filesystem
$ sudo chown -R 105.105 $PWD/cartodb_pgdata

# 这样即使 container 重启，之前建立的数据集也不会丢失了。
$ docker run -d -p 80:80 -h ip -v $PWD/cartodb_pgdata:/var/lib/postgresql sverhoeven/cartodb
```
