# Kraken

## 1. 安装

```bash
$ conda create -n kraken kraken2 bracken=2.2
$ conda activate kraken
(kraken)$ wget https://www.ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz
(kraken)$ tar zxvf minikraken2_v2_8GB.tgz -C $HOME/dbs/minikraken2
```

## 2. 使用

```bash
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report.txt assembly.fna > result.kraken
(kraken)$ kraken2 --use-names --threads 4 --db $HOME/dbs/minikraken2 --report report2.txt --fastq-input --paired S1_R1.fq.gz S1_R2.fq.gz > result2.kraken
(kraken)$ wget https://ccb.jhu.edu/software/bracken/dl/minikraken2_v2/database200mers.kmer_distrib
(kraken)$ mv database200mers.kmer_distrib $HOME/dbs/minikraken2
(kraken)$ bracken -d $HOME/dbs/minikraken2 -i kraken2.report -o bracken.species.txt -l S
```

## 3. 示例
