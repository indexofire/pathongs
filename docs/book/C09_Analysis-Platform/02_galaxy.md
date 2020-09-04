# Galaxy 快速部署方式

---

![banner](../assets/images/9.2/banner.jpg)

## 1. 小型实验室快速解决方案

使用 Amazon EC2 对于不熟悉 Linux 和服务器配置的用户而言是很方便，节省用户在硬件的经济投入，配置学习以及维护能力的培养。但对于测序数据不适合上传到第三方服务器的实验室或要结合本地 LIMS 系统的实验室，要快速构建 Galaxy 以及相配套的其他软件，目前来说最简便的方法之一是利用 docker 部署到本地服务器上。

### 1.1 安装 Docker

Ubuntu 14.04 已经自带了 docker 安装包，不过版本相对较老，想使用 docker 最新版的功能，添加源来安装 docker:

```bash
# 添加 docker 源
$ sudo apt-get install apt-transport-https
$ sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 36A1D7869245C8950F966E92D8576A8BA88D21E9
$ sudo bash -c "echo deb https://get.docker.io/ubuntu docker main > /etc/apt/sources.list.d/docker.list"

# 安装 docker
$ sudo apt-get update
$ sudo apt-get install lxc-docker

# 安装完成后启动 docker 服务
$ sudo service docker start
```

### 1.2 Docker入门

首先要对 docker 镜像和容器的概念有所了解。细节可以参阅 [Docker 从入门到实践](https://yeasy.gitbooks.io/docker_practice/content)

```bash
$ docker run ubuntu:14.04 /bin/echo 'Hello world'
```

### 1.3 安装 galaxy docker 镜像

galaxy 的 docker 镜像可以自己来创建，建议使用`docker galaxy-stable`，源代码可以在 [Github](https://github.com/bgruening/docker-galaxy-stable) 下载。也可以直接到官方 Hun Registry 里下载`galaxy-stable`镜像

```bash
$ sudo docker pull bgruening/galaxy-stable
```

官方 Registry 非常慢，可以通过国内 daocloud 服务商的加速器来加快`docker pull`过程：首先注册一个 [daocloud](daocloud.io) 用户，然后在命令行中添加 registry mirror。

```bash
# ****** 为daocloud分配给你的ID
$ echo "DOCKER_OPTS=\"\$DOCKER_OPTS --registry-mirror=http://******.m.daocloud.io\"" | sudo tee -a /etc/default/docker
$ sudo service docker restart
$ sudo docker pull bgruening/galaxy-stable
```

该脚本可以将`--registry-mirror`加入到你的Docker配置文件`/etc/default/docker`中。适用于Ubuntu14.04，其他版本根据docker配置文件和环境变量，可能有细微不同。

```bash
$ wget https://github.com/bgruening/docker-galaxy-stable/archive/15.10.tar.gz
$ tar zxf 15.10.tar.gz
$ cd docker-galaxy-stable-15.10
```

### 1.4 启动 galax-stable 容器

`docker run -d`表示以daemon方式运行docker，执行该命令后，docker完成galaxy初始化需要花几分钟时间才能访问。如果是本地运行，用浏览器访问`127.0.0.1:8080`可以看到galaxy首页，如果是服务器上用ssh连接执行命令，要访问服务器IP加上8080端口。

```bash
# 后台运行 galaxy 容器
$ sudo docker run -d -p 8080:80 -p 8021:21 bgruening/galaxy-stable

# 交互方式运行 galaxy 容器
$ sudo docker run -i -t -p 8080:80 bgruening/galaxy-stable /bin/bash
root$ sh run.sh
```

有时候需要进入docker容器中进行操作，就可以以`docker run -i`交互模式进行访问。进入容器后运行`sh run.sh`可以DEBUG方式运行galaxy，适合本地测试使用。

### 1.5 复制容器内文件

```bash
# fc3e62e0471d 是想要获取文件的所在容器。foo.txt是想要获得的文件。
$ sudo docker cp fc3ea62e471d:/home/foo.txt .
```

### 1.6 加载数据卷

需要分析的数据通过添加外部数据卷来实现。

```bash
# 添加服务器上的 /mydata 卷到容器中
$ sudo docker create -v /mydata --name my_data_vol bgruening/galaxy-stable /bin/bash

# 或者在运行时将本地卷`/mydata`加入到容器中`/container_data`位置
$ sudo docker run -d -p 8080:80 -v /mydata:/container_data/ bgruening/galaxy-stable
```

镜像是只读的，当`Ctrl+D`方式退出容器后，再次进入容器时你上次以添加的内容是看不到的。如果想要从上一次运行的容易中获得文件可以用`docker cp`的方法，不过你得记住上一次运行的container id号。

### 1.7 删除所有不运行的容器

```bash
$ sudo docker ps -a | cut -d ' ' -f 1 | sudo xargs docker rm
```

### 1.8 删除镜像

```bash
# IMAGE ID 是该镜像的ID，如果镜像还有容器运行，或有其他镜像的依赖关系，则无法删除要先删除容器或其他镜像。
$ sudo docker rmi IMAGE_ID
```

---

## 2. 本地安装与基本配置

本节介绍 Galaxy 的最基本下载安装与使用。最适合的场景为个人电脑，单用户使用的情况。

### 2.1 下载与安装

Galaxy 作为一款开源软件，其代码库托管在 http://bitbucket.org ，先安装 mercurial ，然后用 hg 工具将 galaxy 代码库克隆到本地。

```bash
$ sudo apt-get install mercurial
$ hg clone https://bitbucket.org/galaxy/galaxy-dist/
$ cd galaxy-dist
```

可以用 `hg branch` 命令查看代码分支是否为 stable，如果是其他分支（galaxy代码库有另一分支'default'），切换到 stable 分支：`hg update stable`。在生产环境下建议使用 stable 分支的代码。

### 2.2 配置与运行

克隆到本地的 stable 代码，一般Linux系统自带Python就可以直接运行了。

```bash
$ ./run.sh
```

运行 `run.sh`，这个shell脚本程序会自动完成初始化数据，依赖库下载，数据库迁移等一系列操作，当看到终端显示`serving on http://127.0.0.1:8080`时，可以打开浏览器，访问 http://127.0.0.1:8080 即可看到galaxy的界面。

### 2.3 添加官方 toolshed 中的工具

默认的galaxy只带有基本的工具，对于实际工作中需要的各种分析软件，需要添加到自己建立的 galaxy 实例中。

高通量测序的生物信息学软件大多是基于命令行的开源工具。galaxy利用python语言将这些工具粘合到galaxy实例中，使得用户可以在web界面中直接调用命令行工具对数据进行操作。

galaxy有一个toolshed（工具库）的概念：`https://toolshed.g2.bx.psu.edu/`（官方维护的toolshed），许多著名的工具已经被移植到toolshed中，可以直接被安装到galaxy里，此外也有许多第三方的toolshed包可以添加，甚至掌握了一些python脚本和galaxy xml规范后，也可以自己添加一些分析工具到galaxy中。

除了分析软件外，toolshed还包含创建的数据类型，以及工作流等。

首先修改配置文件 `tool_sheds_conf.xml`

```bash
$ cp config/tool_sheds_conf.xml.sample config/tool_sheds_conf.xml
```

其次在galaxy.ini中配置依赖包的安装目录（上一部分已经添加了这个参数，将依赖包安装在`tool_dep`），然后添加管理员帐号，比如你之前用admin@localhost.com注册的galaxy实例，就将galaxy.ini中`admin_users`设置为：

```bash
admin_users = admin@localhost.com
```

重启你的galaxy实例，用admin@localhost.com用户登陆，你就有权限访问`http://127.0.0.1:8080/admin`，在admin界面可以看到`Tool sheds`下有`Search and browse tool sheds`链接，点击后可以看到默认的2个toolsheds源。

进入`Galaxy main tool shed`，工具列表上访有搜索框，在这里输入你要安装的工具名称比如`spades`后，回车进行检索。

![Instance](../assets/img/appendix_a5_1.png)

点击结果列表中的`spades`下拉菜单，选择`Preview and install`，在转向页面中点击`Install to Galaxy`后，会出现如下图的提示。

![Instance](../assets/img/appendix_a5_2.png)

在`add new tool panel section`中输入`Assembly`，将`spades`工具归类到`Assembly`这个新建的工具类别中，点击页面底部的`Install`按钮开始安装。

---

## 3. 构建单机产环境

作为个人尝试，前面的步骤在PC机上已经可以正常运行使用。对于要作为生产环境下多用户使用，建议使用专门的代理服务器和数据库来增强效率和速度。这里采用nginx+postgresql构建生产环境下的 Galaxy 服务。这里只是简单的介绍一下最基本的配置方式，对于高负载的web server设置又是另一个很复杂的话题了，这里就不具体展开。也可以参考别人做的[galaxy dockerfile](https://registry.hub.docker.com/u/bgruening/galaxy-stable/dockerfile/)。

### 3.1 系统设置

首先在ubuntu下新建一个用户galaxy。

```bash
$ sudo adduser galaxy
```

按照提示，在弹出提示符输入相应内容，主要填好密码即可，其他可以留空。然后切换到galaxy用户：

```bash
$ su galaxy
$ cd
```

重复`Galaxy 本地安装与配置`的第一步`下载与安装`步骤,建立配置文件：

```bash
$ cp config/galaxy.ini.sample config/galaxy.ini
$ vim config/galaxy.ini
```

将配置文件galaxy.ini设置如下，你也可以将内容直接保存成galaxy.ini

```bash
[server:main]
use = egg:Paste#http
host = 0.0.0.0
use_threadpool = True
threadpool_kill_thread_limit = 10800

[filter:gzip]
use = egg:Paste#gzip

[filter:proxy-prefix]
use = egg:PasteDeploy#prefix
prefix = /galaxy

[app:main]
paste.app_factory = galaxy.web.buildapp:app_factory
database_connection = postgresql://galaxy:galaxy@localhost:5432/galaxyserver
tool_dependency_dir = tool_dep
use_nglims = False
nglims_config_file = tool-data/nglims.yaml
debug = False
use_interactive = False
admin_users = admin@localhost.com
```

### 3.2 安装 postgresql

安装postgresql并建立galaxy数据库，这里用户名和密码都设置为`galaxy`，要与 galaxy.ini 中`database_connection`参数对应的值一致。

```bash
$ sudo apt-get install postgresql-9.3
$ su - postgres
$ psql template1

> CREATE USER galaxy WITH PASSWORD 'galaxy';
> CREATE DATABASE galaxyserver;
> GRANT ALL PRIVILEGES ON DATABASE galaxyserver to galaxy;
> \q
```

### 3.3 安装 nginx

用nginx做反向代理，处理请求。

```bash
$ sudo apt-get install nginx
```

设置nginx.conf

```bash
http {
    upstream galaxy_app {
        server localhost:8080;
    }

    server {
        client_max_body_size 10G;
        location / {
            proxy_pass   http://galaxy_app;
            proxy_set_header   X-Forwarded-Host $host;
            proxy_set_header   X-Forwarded-For  $proxy_add_x_forwarded_for;
        }
    }
}
```

## Rreference

1. [Docker 从入门到实践](https://yeasy.gitbooks.io/docker_practice/content)
2. [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable)
