# 下载BioSample信息

!!! Abstract "内容简介"

    获得基因组后，希望了解菌株的meta信息，就需要查看biosample数据库中基因组对应的菌株信息。

从获取一个已知菌株在biosample数据库中的记录开始。

```bash
# 查询菌株ATCC25922的biosample
$ esearch -db biosample -query "ATCC25922"
<ENTREZ_DIRECT>
  <Db>biosample</Db>
  <WebEnv>MCID_661d37307a915242cb1612a2</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>42</Count>
  <Step>1</Step>
</ENTREZ_DIRECT>
```

结果可以获得42个查询结果，将结果信息生成纯文本格式查看。

```bash
# 获得纯文本格式
$ esearch -db biosample -query "ATCC25922" | efetch -format txt
1: OneHealth
Identifiers: BioSample: SAMN40301290; Sample name: WVDA_M07713_Ecoli_ATCC25922; SRA: SRS20722933
Organism: Escherichia coli
Attributes:
    /strain="Escherichia coli_Lot 705489"
    /collected by="West Virginia Department of Agriculture"
    /collection date="2024-02-27"
    /geographic location="USA:WV"
    /isolation source="new culture swab"
    /source type="other"
    /purpose of sampling="research [GENEPIO:0100003]"
    /project name="GenomeTrakr; LFFM-FY4"
    /Interagency Food Safety Analytics Collaboration (IFSAC) category="clinical/research"
    /sequenced by="West Virginia Department of Agriculture"
    /LexMaprStandardizedIsolationSource="new culture swab"
Accession: SAMN40301290	ID: 40301290
...
```

选择我们需要的基因组，并下载biosample信息。

```bash
$ esearch 
```