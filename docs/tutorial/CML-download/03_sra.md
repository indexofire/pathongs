# 下载raw测序数据

!!! Abstract "内容简介"

    本节介绍下载NCBI/EBI/DDBJ的高通量测序原始数据。比如SRA数据，或者以fastq.gz储存的原始数据。


## Kingfisher工具



### 安装 Kingfisher

`aspera-cli`工具是IBM开发的aspera命令行工具，ascli是该工具运行命令，同时软件包内含ascp, curl等工具，因此安装该工具即可在命令行下载aspera服务器端数据。aspera_bypass_dsa.pem

#### conda虚拟环境方式

```shell
# 创建虚拟环境
$ mamba create -n getraw
$ mamba activate getraw
# 安装软件
(getraw)$ mamba install kingfisher

# 采用aspera下载，需要安装aspera-cli
(getraw)$ mamba install aspera-cli
# 通过aws云端数据下载，需要安装
(getraw)$ mamba install awscli
# 通过google-cloud云端数据下载，需要安装
(getraw)$ mamba 
```

!!! tip "注意事项"

    conda环境如果默认添加了bioconda和conda-forge的channels，安装aspera-cli会搜索到v4版本的aspera-cli，IBM用ruby重写了控制界面，可以使用ascli工具来调用下载。当然也可以直接使用ascp程序进行下载。而之前的v3版本需要通过channels：hcc来下载(mamba install -c hcc aspera-cli)。

    v3和v4版本最大的区别在与使用的公钥文件名是不一样的：aspera_bypass_dsa.pem(v4)和asperaweb_id_dsa.openssh(v3)。而kingfisher在软件包内直接提供了openssh文件。

#### docker虚拟环境方式

```shell
# 
$

```

### 使用 Kingfisher 下载

```shell
# 下载SRA登陆号为SRR17840141，下载方式为aspera，从ENA下载
(raw-dl)$ kingfisher get -r SRR17840141 -m ena-ascp

```

### Kingfisher 问题或缺点

1. 当使用`--run-identifiers-list`批量下载时，使用`ena-ascp`模式下载时容易报错STDERR，导致下载中止。因此



---

## nf-core/fetchngs