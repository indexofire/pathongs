# 命令行下载 NCBI Assembly 数据

!!! Abstract "内容简介"
    
    通过命令行下载[Assembly][]数据的方法和工具，不同工具之间的优缺点，可以根据自己网络情况和需求选择不同的下载方式，本节介绍研究工作需求下载数据时常用的一些方法。
    
    本节使用的工具包括：

    - [entrez-direct][]: v16.2
    - [ncbi-datasets-cli][]: v16.11.0
    - [aria2][]: v1.37.0
    - [curl][]: v8.5.0
    - [jq][]: v1.6

## 下载单个基因组数据

下载单独的[Assembly][]数据很简单，一般通过网页端即可。在终端下载可以使用的工具和方式：

**FTP下载**

如果已知assembly数据库的某个基因组accession，比如GCF_000006745.1，可以检索获得ftp url信息，通过wget,curl或aria2c等工具进行下载。

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
$ datasets download genome accession GCF_000006745.1
$ 7z x ncbi_dataset.zip
$ ls ncbi-dataset/data/GCF_000006745.1
GCF_000006745.1_ASM674v1_genomic.fna

# 下载其他格式数据，比如下载gbk,gff格式的数据
$ datasets download genome accession --include gbff,gff3 GCF_000006745.1
$ 7z x ncbi_dataset.zip
$ ls ncbi-dataset/data/GCF_000006745.1
GCF_000006745.1_ASM674v1_genomic.gbff
GCF_000006745.1_ASM674v1_genomic.gff
```

下载某个菌种的参考基因组

```bash
# 下载霍乱弧菌参考基因组
$ datasets download genome taxon 666 --reference
```

## 下载多个基因组数据

**通过API搜索需要的数据**

通过WEB访问[NCBI][]的网页，利用其[Entrez][]工具检索到所需要的基因组数据，通过http方式下载是普通用户的一般方式。当在服务器终端下，我们可以利用[NCBI][]提供的[Entrez][entrez-direct]命令行工具检索并下载数据。

**以检索条件下载**

[entrez-direct][]现在命令行根据筛选条件进行检索，获得ftp地址后通过[aria2][]工具下载。

```bash
# 检索2023年发布的refseq基因组数据
$ esearch -db assembly -query "("salmonella"[Organism:exp] AND ("2023/01/01"[SeqReleaseDate]:"2023/12/31"[SeqReleaseDate])) AND "latest+refseq"[filter]"
<ENTREZ_DIRECT>
  <Db>assembly</Db>
  <WebEnv>MCID_661cbbe57c9c0a08887ae7ca</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>1174</Count>
  <Step>1</Step>

# 下载1174个命中的结果信息
$ esearch -db assembly -query "("salmonella"[Organism:exp] AND \
> ("2023/01/01"[SeqReleaseDate]:"2023/12/31"[SeqReleaseDate])) \
> AND "latest+refseq"[filter]" | \
> efetch -format docsum > result.xml
$ xtract -input result.xml -pattern DocumentSummary -element FtpPath_RefSeq | \
> awk -F’/’ ‘{print $0”/”$NF”_genomic.fna.gz”}’ > ftp_url
$ aria2c -j 3 -i ftp_url
```

通过检索得到的xml信息中accession信息生成列表清单，再通过datasets进行下载。dehydrated模式可以获得多个基因组下载地址，并通过rehydrated方式快速批量下载，效果好于ftp下载方式。是目前大规模下载[Assembly][]基因组数据最好的方法。

```bash
# 提取之前1174个命中的基因组的accession
$ xtract -input result.xml -pattern DocumentSummary \
> -element AssemblyAccession > accession
# 生成下载信息，包含在压缩包ncbi_datasets.zip
$ datasets download genome accession --inputfile accession --dehydrated
# 解压缩文件，下载信息保存在fetch.txt中
$ 7z x ncbi_datasets.zip
# rehydrated，下载基因组
$ datasets rehydrated --directory .
```

我们还可以使用`datasets summary`命令进一步筛选基因组。

```bash
# 使用summary命令生成json格式的基因组信息
$ datasets summary genome accession --inputfile accession > result.json
# jq工具进行筛选，保留checkM值大于99.7的基因组
$ cat result.json | jq -r '.[][] | .accession, .checkm_info.completeness' | \
> xargs -n2 | awk '{if($2>99.7)print}' > filter
$ wc -l filter
199

# 剩下符合条件基因组199个，dehydrated方式下载
$ datasets download genome accession --inputfile filter --dehydrated
$ 7z x ncbi_datasets.zip
$ datasets rehydrated --directory .
```

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

采用datasets下载物种[Assembly][]数据

```bash
$ datasets download genome taxon $txid --dehydrated
$ 7z x ncbi_datasets.zip
$ datasets rehydrated --directory .
```

!!! tip

    NCBI最近即将更新其Genome和Assembly页面，新的版本和原页面功能有较大改变。

## 第三方工具

其他常用的第三方小工具：

- [get_assemblies](https://github.com/davised/get_assemblies)
- [ncbi-kit](https://github.com/jameslz/ncbi-kit)

个人使用感受，目前下载[Assembly][]数据库的基因组数据最优秀的工具还是官方的[ncbi-datasets-cli][]工具，起`rehydrated`下载方式高效且不易出错。

[Assembly]: https://www.ncbi.nlm.nih.gov/assembly/
[entrez-direct]: https://www.ncbi.nlm.nih.gov/books/NBK179288/
[ncbi-datasets-cli]: https://www.ncbi.nlm.nih.gov/datasets/
[aria2]: https://aria2.github.io/
[curl]: https://curl.se
[jq]: https://jqlang.github.io/jq/
[NCBI]: https://www.ncbi.nlm.nih.gov/
[Entrez]: https://www.ncbi.nlm.nih.gov/search/