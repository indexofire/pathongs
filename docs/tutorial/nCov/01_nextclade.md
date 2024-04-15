# NextClade 分析病毒变异位点

!!! note "内容简介"

    nextclade是用于新冠等病毒突变位点溯源的工具，分为命令行和在线2个版本。在线版本如clades.nextstrain.org是一个实例。如果针对大规模数据，可以直接使用命令行版本来进行分析。

## nextclade-cli

**安装**

```shell
# 安装nextclade
(sars)$ mamba install nextclade=3.3.1
```

**使用**

```shell
# 下载参考数据库
(sars)$ nextclade dataset get \
> --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' \
> --output-dir 'data/sars-cov-2'
```

已完成测序数据: run/*.fasta

```shell
# 基因组数据分析
(sars)$ nextclade run \
> --input-dataset data/sars-cov-2 \
> --output-all=output/ \
> run/*.fasta

# 查看生成数据
(sars)$ ls output
nextclade.aligned.fasta
nextclade.auspice.json
nextclade.csv
nextclade.cds_translation.ORF3a.fasta
nextclade.cds_translation.ORF6.fasta
nextclade.cds_translation.E.fasta      
nextclade.cds_translation.ORF7a.fasta  
nextclade.cds_translation.M.fasta      
nextclade.cds_translation.ORF7b.fasta  
nextclade.cds_translation.N.fasta      
nextclade.cds_translation.ORF8.fasta   
nextclade.cds_translation.ORF1a.fasta  
nextclade.cds_translation.ORF9b.fasta
nextclade.cds_translation.ORF1b.fasta  
nextclade.cds_translation.S.fasta
nextclade.json
nextclade.ndjson
nextclade.nwk
nextclade.tsv

# 可以只生成单独的数据格式，比如只生成tsv格式。
(sars)$ nextclade run \
> --input-dataset data/sars-cov-2 \
> -t output/output.tsv \
> run/*.fasta

# 查看具体结果
(sars)$ less output/nextclade.tsv
```

## nextclade-web

对于国内疾控实验室来说，涉及到敏感数据，往往不建议之间上传数据到公共服务器上，因此在本地部署nextclade-web，提供可视化界面给研究人员进行毒株基因组分析提供了很好的解决方案。

**安装**

```shell
# 克隆nextclade代码仓库
(sars)$ git clone https://github.com/nextstrain/nextclade
(sars)$ cd nextclade

# 安装依赖工具
# 1. 安装rust
(sars)$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# 按照提示完成安装后，将路径添加到环境变量中
(sars)$ source ~/.bashrc
# 2. 安装wasm
(sars)$ cargo install wasm-pack
# 3. 安装nodejs
(sars)$ mamba install nodejs
# 4.安装yarn
(sars)$ npm install -g yarn

# 安装nextclade-web
(sars)$ cd packages/nextclade-web
(sars)$ yarn install
(sars)$ yarn wasm-prod

# 运行测序环境
(sars)$ yarn dev

# 打开浏览器，访问localhost:3000
(sars)$ firefox localhost:3000
```

**使用**

nextclade-web基本使用方法与在线工具[nextclades](clades.nextstrain.org)类似。




通过构建自定义数据库，可以将区域毒株联合分析，不仅可以构建地区数据库，而且一旦发生疫情可以快速溯源到近缘毒株以及所有SNVs位点，便于流行病学溯源。

```shell
(sars)$ 
```



[参考文档](https://github.com/nextstrain/nextclade/blob/master/docs/dev/developer-guide.md)