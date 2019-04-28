# 二代+三代测序数据混合拼接

## Spades

```bash
# Spades 混合拼接同一个样本 miseq 和 nanopore 测序数据
$ spades.py -m -k 33,55,77,99,121 --only_assembler -1 S1_R1.fastq.gz -2 S1_R2.fastq.gz --nanopore minion.fastq.gz -o hybrid -t 40
```

## Unicycle

```bash
# 去除接头，将fastq文件压缩
(nanopore)$ for fq in {0..11}; do porechop -t 8 -i *_${fq}.fastq -o ${fq}.fastq.gz; done

# 去除低质量碱基，将reads平均Q值低于9的去除
(nanopore)$ for trim in {0..11}; do zcat ${trim}.fastq.gz | NanoFilt -q 9 >> trimmed.fastq; done

# unicycler 混合组装illumina和nanopore数据
(nanopore)$ unicycler -1 R1.fastq.gz -2 R2.fastq.gz -l trimmed.fastq.gz -o hybrid
```
