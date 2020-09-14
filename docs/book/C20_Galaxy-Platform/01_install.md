# 安装 Galaxy

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 安装方式1.

下载源代码

```bash
# git clone github上托管的源代码
# -b 指定分支，比如目前最新版本的分支release_20.05，追求stable稳定发布版的可以选择这个分支
$ git clone -b release_20.05 https://github.com/galaxyproject/galaxy/

# 如果尝试新功能的，可以下载master分支代码
$ git clone -b master https://github.com/galaxyproject/galaxy/
```

运行程序，最理想的方式还是采用虚拟环境，比如用conda建立。因为galaxy是python开发的，也可以利用众多的python虚拟环境包建立。

```bash
# conda 环境
$ conda create -n galaxy
$ conda activate galaxy

# 初始化galaxy，第一次运行run.sh脚本会检查环境并安装依赖软件
(galaxy)$ cd galaxy && sh run.sh
```

打开浏览器，访问`127.0.0.1:8080`可以看到galaxy的显示页面。

## 安装方式2. docker方式

Galaxy 由于依赖众多，安装比较复杂，推荐用 docker 镜像安装，比较简单。但是对docker需要有一定的了解，特别是配置数据卷方面。

```bash
# 下载镜像
$ docker pull bgruening/galaxy-stable

# 运行images
$ docker run -d -p 8080:80 -p 8021:21 -p 8022:22 bgruening/galaxy-stable
# 打开浏览器，访问 http://127.0.0.1:8080 可以看到 galaxy 界面
```





## Reference

1. https://mimodd.readthedocs.io/en/latest/install_galaxy.html
2.
