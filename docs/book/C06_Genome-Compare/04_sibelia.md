# C-Sibelia

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

[Sibelia][] 是一个由 [St. Petersburg Academic University of the Russian Academy of Sciences][] 开发的主要用于病原微生物的比较基因组学工具。它可以用来发现微生物的共有区域，非共有区域以及重组区域等。通过对共线性区域的分析与比较，并采用 [circos][] 以及 [d3js][] 生成可视化结果。

[Sibelia][] 有2个组件：一个叫`Sibelia`，可以用来处理多个基因组完成图的共线性分析。另一个叫`C-Sibelia`，可以用来比较2个基因组（基因组scaffolds与参考基因组）。

## 1. 安装 Sibelia

### 1.1 安装依赖包

如果要自己编译安装[Sibelia][]，需要先安装以下的 C++ 库和软件。然后自行编译。

```bash
# install libdivsufsort
$ wget https://github.com/y-256/libdivsufsort/archive/2.0.1.tar.gz
$ tar zxf 2.0.1.tar.gz
$ cd libdivsufsort-2.0.1
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" ..
$ make
$ sudo make install

# install tclap
$ wget http://sourceforge.net/projects/tclap/files/tclap-1.2.1.tar.gz
$ tar zxf tclap-1.2.1.tar.gz
$ cd tclap-1.2.1
$ ./configure
$ make
$ sudo make install

# install seqan
$ sudo apt-get install lib
$ wget
http://packages.seqan.de/seqan-library/seqan-library-2.0.1.tar.bz2
$ tar xjf seqan-library-2.0.1.tar.bz2
$ sudo cp -R seqan-library-2.0.1/include/seqan /usr/include
$ sudo cp -R seqan-library-2.0.1/share/doc/seqan /usr/share/doc

# install boost
$ sudo apt-get install libboost-dev

# build Sibelia
$ git clone https://github.com/bioinf/Sibelia
$ cd Sibelia/build
$ cmake ../src
$ make
$ sudo make install
```

### 1.2 下载 Sibelia 预编译包

```bash
$ wget http://sourceforge.net/projects/sibelia-bio/files/3.0.6/Sibelia-3.0.6-Linux.tar.gz
$ tar zxf Sibelia-3.0.6-Linux.tar.gz -C ~/app
$ echo 'export PATH="$PATH:$HOME/app/Sibelia-3.0.6/bin"' >> ~/.bashc
```

## 2. 使用 C-Sibelia

### 2.1 比较测序结果与参比基因组

```bash
# run sibelia
$ C-Sibelia.py reference.fasta sequencing.fasta -o output

# plot circos image
$ circos -conf output/circos/circos.conf
# view circos result
$ evince output/circos/circos.svg
# or view the lines plot
$ google-chrome output/d3_blocks_diagram.html
```

## 3. Reference

1. [Sibelia](http://bioinf.spbau.ru/sibelia)


[Sibelia]: http://bioinf.spbau.ru/sibelia
[St. Petersburg Academic University of the Russian Academy of Sciences]: http://spbau.ru
[circos]: http://circos.ca
[d3js]: http://d3js.org
