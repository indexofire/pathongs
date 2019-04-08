# GI to Accession Number



NCBI 已经取消了`gi`，全面转向accession。如果手上有数据还是`gi`号，建议将其进行转换。

## 1. 转换方法

### 1.1 使用E-utilies

打开浏览器，地址栏里输入以下内容，将你想查询`gi`编号放在`id=`后面，多个`gi`可以用逗号分隔。访问该地址即可获得对应的 accession。
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=663070995,568815587&rettype=acc

### 1.2 使用edirect

```bash
# efetch 获取 accession
$ efetch -db nuccore -id 663070995 -format acc

# 批量修改 VFDB 文件，无NCBI账户每秒只能发出3次申请
$ cp VFDB_setB_nt.fas VFDB_setB_nt.fas.bak
$ for i in $(grep 'gi:' VFDB_setB_nt.fas | awk -F'(' '{print $2}' | awk -F')' '{print $1}'); do r=`efetch -db nuccore -id ${i#gi:} -format acc`; sed -i 's/'$i'/'$r'/g'; sleep 0.1; done
```

### 1.3 gi2accession.py

下载 `ftp://ftp.ncbi.nlm.nih.gov/genbank/livelists/gi2acc_mapping/` 下数据库和修改脚本gi2_accession.py。数据库文件比较大，除非需要转换的量大或者频率高，一般不建议

```bash
$ python gi2_accession.py
gi: 42
42  CAA44840.1  416

# 批量获得结果
$ python gi2_accession.py < list_ids.txt
```

## 2. 补充材料

### 2.1 Refseq Accession 含义

| Accession prefix   | Molecule type |   Comment   |
| :----- | :-------: | :------ |
| AC_    |   Genomic |  Complete genomic molecule, usually alternate assembly  |
| NC_    |   Genomic |  Complete genomic molecule, usually alternate assembly  |
| NG_    |   Genomic |  Incomplete genomic region  |
| NT_    |   Genomic |  Contig or scaffold, clone-based or WGS*  |
| NW_    |   Genomic |  Contig or scaffold, primarily WGS[1]  |
| NZ_[2]  |   Genomic |  Complete genomes and unfinished WGS data  |
| NM_    |   mRNA    |  Protein-coding transcripts (usually curated)  |
| NR_    |   RNA     |  Non-protein-coding transcripts  |
| XM_[3] |   mRNA    |  Predicted model protein-coding transcript  |
| XR_[3] |   RNA     |  Predicted model non-protein-coding transcript  |
| AP_    |   Protein |  Protein	Annotated on AC_ alternate assembly  |
| NP_    |   Protein |  Protein	Associated with an NM_ or NC_ accession  |
| YP_[3] |   Protein |  Annotated on genomic molecules without an instantiated transcript record  |
| XP_    |   Protein |  Predicted model, associated with an XM_ accession  |
| WP_    |   Protein |  Non-redundant across multiple strains and species  |
- [1]：whole genomic sequence
- [2]：An ordered collection of WGS sequence for a genome.
- [3]：


## Reference

1. [Refseq Accession Dict](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly)
2. [NCBI Youtube Tutorial](https://www.youtube.com/watch?v=rIDQEnnOr6g)
