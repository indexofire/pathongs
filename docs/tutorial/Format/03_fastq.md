# FastQ 格式

---

fastq 格式的数据由若干条reads组成，每条read由4行表示信息为

第一行内容包括测序仪信息，实验信息和序列信息等

```
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
```




第二行内容是序列

```

```

第三行为标识符号，没有特殊含义

```

```

第四行为对应的序列质量

```

```
