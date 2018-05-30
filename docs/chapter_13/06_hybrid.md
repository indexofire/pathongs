# Hybrid assemble

## Spades

```bash
# Spades 混合拼接同一个样本 miseq 和 nanopore 测序数据
$ spades.py -m -k 33,55,77,99,121 --only_assembler -1 S1_R1.fastq.gz -2 S1_R2.fastq.gz --nanopore minion.fastq.gz -o hybrid -t 40
```
