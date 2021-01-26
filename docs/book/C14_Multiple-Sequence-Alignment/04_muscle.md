# Muscle

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abtract "内容简介"
    Muscle 是一个老牌的多重序列比对工具，它是一个progressive alignment工具，精度与mafft相似，但速度更快。[参考文献](https://academic.oup.com/nar/article/32/5/1792/2380623)

Muscle 构建多重序列比对基本步骤是3部分：草图Progressive式比对（利用kmer进行距离比较）；精细Progressive比对（Kimura distance重评估）；精细修正

## 安装

使用版本：

- muscle: v3.8.1151

```bash
# muscle 没有需要其他依赖，可以直接安装在默认的环境中
$ conda activate seqs
(msa)$ conda install muscle
```

## 使用

```bash
# seqs.fas 是 multifasta 文件，seqs.aln 是多重序列比对结果
(msa)$ muscle -in seqs.fas -out seqs.aln

# 对于近源序列，可以添加 -diags 参数，并减少回调次数加快运行
(msa)$ muscle -in seqs.fas -out seqs.aln -diags -maxiters 2
```

## Reference

1. [User Guide](http://www.drive5.com/muscle/muscle_userguide3.8.html)
