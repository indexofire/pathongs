# Roary

## 使用

**输入数据**

Roary 的输入端需要 gff 格式的数据文件。由于 NCBI 下载的 gff 是不含核酸序列的 gff 格式文件，无法直接用来分析。Prokka 生成的 gff 格式文件包含核酸序列，所以可以下载 NCBI 上的 fna 文件然后用 prokka 注释后再用 roary 分析。

```bash
# 修改路径
$ sed -i 's/.gbff/.fna/g' ftp-path.txt
$ wget -i ftp-path.txt --no-passive-ftp
$ gunzip *.gz

# 将 Assembly 文件名修改为菌株名称
$ for i in *.fna; do mv $i `head -1 $i | awk '{print $5, $6}' | awk -F',' '{print $1".fna"}' | sed 's/CMCC /CMCC-/g'`; done
# prokka 进行注释并以菌株名为文件夹名输出
$ for i in *.fna; do k=`echo $i | sed 's/.fna//g'` && prokka $i --prefix $k -outdir $k --cpus 40; done
# 复制出 gff 文件到 gff_folder 中，由 roary 进行分析
$ cp */*.gff gff_folder/; done
```

**使用方法**

roary 用 blastp 对 gff 文件中的序列进行 orthologs 分析，获得 pangenome 和 coregenome 结果，还可以使用 mafft 对核心基因组进行序列比对，生成系统发生树。

```bash
$ roary -e --mafft -p 40 gff_folder/*.gff -f roary
```

**结果文件**

roary生成的结果文件：

| 结果文件 | 数据内容 |
| -------- | -------- |
| accessory_binary_genes.fa | |
| accessory_binary_genes.fa.newick | |
| accessory_graph.dot | |
| accessory.header.embl | |
| accessory.tab | |
| blast_identity_frequency.Rtab | |
| clustered_proteins | |
| core_accessory_graph.dot | |
| core_accessory.header.embl | embl 格式的文件显示各 accessory 基因|
| core_accessory.tab | accessory 基因在所在的基因组 |
| core_alignment_header.embl | |
| core_gene_alignment.aln | |
| core_gene_alignment.aln.reduced | |
| gene_presence_absence.csv | csv 格式的基因在各个基因组中是否存在的数据文件 |
| gene_presence_absence.Rtab | Rtab 格式的基因在各个基因组中是否存在的数据文件 |
| number_of_conserved_genes.Rtab | Rtab 格式的不同数量基因组所共有基因数 |
| number_of_genes_in_pan_genome.Rtab | Rtab 格式的不同数量基因组的所有基因数 |
| number_of_new_genes.Rtab | Rtab 格式的不同数量基因组所新增的基因数 |
| number_of_unique_genes.Rtab | Rtab 格式的不同数量基因组所特有基因数 |
| pan_genome_reference.fa | |
| summary_statistics.txt | pangenome 分析各种基因数量结果 |
