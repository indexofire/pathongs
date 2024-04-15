# 命令行下载 NCBI Assembly 数据

!!! abstract "内容简介"
    
    命令下有许多下载Assembly数据的方法和工具，根据不同下载方式介绍一些常用的方法。

本教材使用的工具包括：

- entrez-direct
- ncbi-datasets-cli
- aria2c
- curl




## 下载单个基因组数据

**FTP下载**

已知assembly数据库的某个基因组accession，比如GCF_000006745.1，可以检索获得ftp url信息，通过wget,curl或aria2c等工具进行下载。

```bash
# 直接下载GCF_000006745.1基因组fasta数据
$ esearch -db assembly -query 'GCF_000006745.1' | efetch -format docsum | \
> xtract -pattern DocumentSummary -element FtpPath_RefSeq | \
> awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | \
> xargs aria2c
```

**datasets下载**

ncbi-datasets-cli是NCBI推出的下载assembly基因组数据的命令行工具。

```bash
# 默认下载GCF_000006745.1基因组fasta数据
$ datasets download genome accession --filename GCF_000006745.1.fna.gz GCF_000006745.1

```



## 下载多个基因组数据



**通过API搜索**

通过WEB访问NCBI的网页，利用其Entrez工具检索到所需要的基因组数据，通过http方式下载是用户的方式。同样的，当在服务器端，我们可以利用NCBI提供的Entrez命令行工具检索并下载数据。

**以物种为单位进行下载**

需要开展某个物种的研究，希望得到该物种的全部或部分测序数据时，可以从NCBI的taxonomy入手检索会非常方便。因为taxonomy数据库提供种以下的taxon id，这样我们甚至可以对亚种或者型为单位进行检索。以下例子为下载特定血清型沙门菌（乙型副伤寒沙门菌）的Assembly数据库。

```bash
# 检索乙型副伤寒的 taxonomy
$ esearch -db taxonomy -query "salmonella paratyphi B" | efetch -format docsum > taxon.xml
$ txid=$(xtract -input taxon.xml -pattern DocumentSummary -element Id)
# 乙型副伤寒沙门菌的txid
$ echo $txid
57045
```

通过entrez-direct工具获取txid后，可以抓取assembly数据库该txid对应的所有数据。然后将其中基因组在ftp服务器上的地址提取后，采用命令行下载工具如aria2c等进行下载。

!!! tip

    aria2c优点在与支持断点续传，大批量下载的过程发生断联后，未下载完成的文件可以看到aria2c的临时文件，因此可以方便的重新运行aria2c命令以及对应数据的url地址续传数据。

```bash
# 生成乙型副伤寒 Assembly 数据 xml 文件
$ esearch -db assembly -query "txid${id}[Organism:exp] AND latest[filter]" | efetch -format docsum > info.xml
# 提取 Genbank 下载链接
$ xtract -input info.xml -pattern DocumentSummary -element FtpPath_GenBank > ftp.list
# 将info.xml数据中对应基因组的ftp地址提取到文件ftp_url中
$ awk -F’/’ ‘{print $0”/”$NF”_genomic.fna.gz”}’ ftp.list > ftp_url
# aria2c工具下载乙型副伤寒沙门菌的 Assembly 数据，并发连接数3
$ aria2c -j 3 -i ftp_url
# 解压缩文件并重命名
$ gunzip *.gz
# 将文件名重命名为assembly数据库中基因组的accession名称
$ for i in *.fna; do mv $i ${i:0:15}.fna; done
```

!!! tip

    NCBI最近即将更新其Genome和Assembly页面，新的版本和原页面功能有较大改变。

## datasets工具


## ascp工具