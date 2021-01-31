# gggenes 绘制基因分布图

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

```bash
$ prokka genome.fasta --outdir output
$ awk '{print $3 $4 $5}' output/*.tsv > data.tsv 
```

```r
> example_genes$direction <- ifelse(example_genes$strand == "forward", 1, -1)
> ggplot(example_genes, aes(xmin=start, xmax=end, y=molecule, fill=gene, forward=direction)) + 
  geom_gene_arrow() + facet_wrap(~molecule, scale="free", ncol=1) + scale_fill_brewer(palette="Set3") + theme_genes()
```
