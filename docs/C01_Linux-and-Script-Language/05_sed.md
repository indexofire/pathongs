# sed

相比较与awk，sed常用于以行为单位进行替换。

```bash
# 把pattern1替换为pattern2
# 首先搜索pattern1，然后将其替换为pattern2
$ sed s/pattern1/pattern2/g
```
