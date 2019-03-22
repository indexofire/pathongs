# 细菌基因组短序列拼接拼接

微生物基因组测序后的拼接工作，有许多软件可以完成。近几年不断有新的优秀拼接软件涌现，本节仅挑选几个微生物基因组拼接使用较多也是评测中比较优秀的几个工具。也可以先对自己的测序物种做一个 assembly evaluation，选择比较适合的拼接软件。

**常用拼接软件工具列表**

- `SPAdes <http://spades.bioinf.spbau.ru/>`_
- `MaSuRCA <http://www.genome.umd.edu/masurca.html>`_
- Velvet
- A5_miseq

更多工具可以参见`这里 <https://omictools.com/genome-assembly-category>`_


## 1. SPAdes

SPAdes 是由

### 1.1 下载并安装 SPAdes

**下载编译包安装**

下载预先编译的版本包，将可执行文件添加到系统路径中，可以直接调用运行。

```bash
$ curl -O http://spades.bioinf.spbau.ru/release3.11.1/SPAdes-3.9.0-Linux.tar.gz
$ tar zxf SPAdes-3.9.0-Linux.tar.gz -C ~/app
$ cd ~/app
$ sudo ln -s `pwd`/SPAdes-3.11.0-Linux/bin/* /usr/local/sbin
```

**通过conda安装**

```bash
$ conda create -n assembly
$ conda activate assembly
(assembly)$ conda install spades
```

**Docker方式安装**

```bash
$
```

### 1.2 拼装细菌基因组

```bash
# -t 40 表示用40个 threads 进行拼接，根据自己电脑的CPU情况自行设置。
$ spades.py --careful -k 21,33,55,77,99,121 -1 SRR95386_1.fastq -2 SRR955386_2.fastq -o spades_output -t 40
```

SPAdes会尝试不同的Kmer，因此拼装时间也会根据Kmer选择数量成倍增加。对于常见的 Miseq v2/v3 试剂盒，采用 PE150/PE250/PE300 的读长测序，常用的拼接命令是：

```bash
# PE150 读长测序数据
$ spades.py -k 21,33,55,77 --careful -1 SRR95386_1.fastq -2 SRR955386_2.fastq -o SRR95386_output -t 40
# PE250/300 读长测序数据
$ spades.py -k 21,33,55,77,99,127 --careful -1 SRR95386_1.fastq -2 SRR955386_2.fastq -o SRR95386_output -t 40
```

2. Macursa
----------

2.1 下载并安装 MaSuRCA
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   ~/tmp$ curl -O ftp://ftp.genome.umd.edu/pub/MaSuRCA/MaSuRCA-3.1.3.tar.gz
   ~/tmp$ tar zxvf MaSuRCA-3.1.3.tar.gz -C ~/app
   ~/tmp$ cd ~/app
   ~/app$ ./MaSuRCA-2.3.0/install.sh
   ~/app$ sudo ln -s `pwd`/MaSuRCA-2.3.0/bin/masurca /usr/local/sbin

2.2 建立拼装配置文件
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   ~/data$ touch SRR955386_masurca.config
   ~/data$ nano SRR955386_masurca.config

文件内容修改如下：

.. code-block:: bash

   DATA
   PE= p1 500 50 SRR955386_1.fastq SRR955386_2.fastq
   END

   PARAMETERS
   GRAPH_KMER_SIZE=auto
   USE_LINKING_MATES=1
   LIMIT_JUMP_COVERAGE = 60
   KMER_COUNT_THRESHOLD = 1
   NUM_THREADS= 2
   JF_SIZE=100000000
   DO_HOMOPOLYMER_TRIM=0
   END

设置 **GRAPH_KMER_SIZE=auto**，软件会调用Kmer=31来进行拼装。对于MiSeq PE250以上的插片，可以考虑手动设置使用更大的Kmer。

2.3 开始拼装
^^^^^^^^^^^^

.. code-block:: bash

   ~/app$ masurca SRR955386_masurca.config
   ~/app$ ./assemble.sh

--------------------------------------------------------------------------------

3. Velvet
---------

3.1 下载并安装
^^^^^^^^^^^^^^

Velvet 是一个老牌的基因组测序数据拼接软件。Velvet最新版本是1.2.10，可以访问 `官方网站 <https://www.ebi.ac.uk/~zerbino/velvet/>`_ 下载源代码包，也可以通过克隆 `软件仓库 <https://github.com/dzerbino/velvet.git>`_ 的方式获得最新的源代码。

**下载源代码**

普通用户可以下载源代码包自行编译获得软件。

.. code-block:: bash

   ~$ cd /tmp
   ~/tmp$ wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
   ~/tmp$ tar xvf velvet_1.2.10.tgz -C ~/app/velvet

**克隆代码仓库**

如果在软件运行中遇到问题，想试用最新版代码，或是有能力提交issues，或者想改进软件参与开源代码编写的可以选择克隆代码库的方式。

.. code-block:: bash

   ~$ cd ~/tmp
   ~/tmp$ git clone https://github.com/dzerbino/velvet.git

**编译安装**

make可以进行编译，有几个编译参数可以选择，分别是MAXKMERLENGTH, DIRECTORY

.. code-block:: bash

   ~$ cd ~/app/velvet
   ~/app/velvet$ make 'MAXKMERLENGTH=127'
   ~/app/velvet$ sudo cp velveth velvetg /usr/local/sbin

3.2 拼接基因组
^^^^^^^^^^^^^^

velvet软件由2个可执行文件 `velveth` 和 `velvetg` 组成。

.. code-block

   ~$ velvet

3.3 其他组件
^^^^^^^^^^^^

--------------------------------------------------------------------------------

4. A5-miseq
-----------

`A5-miseq <https://sourceforge.net/projects/ngopt/>`__ 是一个用 perl
开发的针对细菌基因组 de novo assembly 的 pipeline
工具。它本身不参与组装，而是通过组合一套工具来完成工作，工具集包括：

- bwa
- samtools
- SGA
- bowtie
- Trimmomatic
- IDBA-UD
- SSPACE

这些工具都以及集成在 A5-miseq 中，不需要另外安装。为了避免不同版本的
samtools，bowtie 对结果产生的差异，建议采用虚拟环境如 docker
等来隔离运行环境。

4.1 安装 A5-miseq
^^^^^^^^^^^^^^^^^

下载预编译包安装
~~~~~~~~~~~~~~~~

.. code-block:: bash

   ~/tmp$ wget http://downloads.sourceforge.net/project/ngopt/a5_miseq_linux_20150522.tar.gz
   ~/tmp$ tar zxvf a5_miseq_linux_20150522.tar.gz -C ~/app
   ~/tmp$ sudo ln -s ~/app/a5_miseq_linux_20150522/bin/a5_pipeline.pl /usr/local/sbin

建立 Docker 容器安装
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   ~/docker$ mkdir -p a5-miseq && cd a5-miseq
   ~/docker/a5-miseq$ touch Dockerfile

.. code-block:: bash

   FROM ubuntu:latest
   MAINTAINER Mark Renton <indexofire@gmail.com>

   RUN apt-get update -qy
   RUN apt-get install -qy openjdk-7-jre-headless file
   ADD http://downloads.sourceforge.net/project/ngopt/a5_miseq_linux_20150522.tar.gz /tmp/a5_miseq.tar.gz
   RUN mkdir /tmp/a5_miseq
   RUN tar xzf /tmp/a5_miseq.tar.gz --directory /tmp/a5_miseq --strip-components=1
   ADD run /usr/local/bin/
   ADD Procfile /
   ENTRYPOINT ["/usr/local/bin/run"]

4.2 使用 A5-miseq
^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ perl a5_pipeline.pl SRR955386_1.fastq SRR955386_2.fastq ~/data/a5_output

4.3 A5-miseq 文档
^^^^^^^^^^^^^^^^^

查看 A5-miseq 工具的使用文档可以用 a5_pipeline.pl 工具查看。

.. code-block:: bash

    # Usage:
    $ a5_pipeline.pl [--begin=1-5] [--end=1-5] [--preprocessed] <lib_file> <out_base>

Reference:
----------

* http://www.chenlianfu.com/?p=1635
* http://www.bbioo.com/lifesciences/40-116878-1.html
* http://www.plob.org/2012/11/21/4797.html
* http://ged.msu.edu/angus/tutorials-2011/short-read-assembly-velvet.html
