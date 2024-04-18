# NextStran平台分析新冠基因组分子流行病学

!!! note "内容简介"

    nextstrain是一个病毒基因组分子流行病学分析平台，可以对毒株基因组的进化关系、流行趋势、地理信息分布、突变位点频度等进行快速分析。

nextstrain工具分为2个核心工具：

- Augur: 生物信息学分析
- Auspice: 数据可视化


## 安装NextStrain

采用conda构建运行环境的安装方式。nextstrain-cli是安装在用户级别，通过创建conda虚拟环境安装平台其他依赖包。

```shell
# 安装nextstrain-cli
$ curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash
$ printf '\n%s\n' 'eval "$("/home/mark/.nextstrain/cli-standalone/nextstrain" init-shell bash)"' >> ~/.bashrc
$ eval "$("/home/mark/.nextstrain/cli-standalone/nextstrain" init-shell bash)"

# 设置conda为默认运行环境
# nextstrain会下载安装micromamba到~/.nextrain路径中
$ nextstrain setup --set-default conda

# 进入nextstrain环境
$ nextstrain shell .
```

## 使用NextStrain分析

**编译ncov新冠数据**

nextstrain是一个开放式的流行病学分析平台，不光可以进行新冠数据分析，还可以针对其他病毒基因组数据。这里以新冠病毒为例，演示

```shell
# 克隆官方代码
$ git clone https://github.com/nextstrain/ncov

```


参见：[官方文档](https://docs.nextstrain.org/projects/ncov)
