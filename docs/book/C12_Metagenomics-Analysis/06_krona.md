# Krona

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

## 安装

```bash
$ conda create -n krona krona
$ conda activate krona
(krona)$ rm -rf $HOME/.conda/envs/krona/opt/krona/taxonomy
(krona)$ mkdir $HOME/dbs/krona/taxonomy
(krona)$ ln -s $HOME/dbs/krona/taxonomy $HOME/.cona/envs/krona/opt/krona/taxonomy
(krona)$ cd $HOME/dbs/krona/taxonomy
(krona)$ ktUpdateTaxonomy.sh .
```
