# Beast2 使用

## 安装

```bash
# 创建 beast2 虚拟环境，安装 beast2
$ conda create -n beast2
$ conda activate beast2
(beast2)$ conda install beast2

# 向conda环境添加变量，优化 java 对 beauti 程序字体渲染
(beast2)$ echo "export _JAVA_OPTIONS='-Dawt.useSystemAAFontSettings=lcd'" >> $CONDA_PREFIX/etc/conda/activate.d/java_home.sh
(beast2)$ echo "unset _JAVA_OPTIONS" >> $CONDA_PREFIX/etc/conda/deactivate.d/java_home.sh
```

