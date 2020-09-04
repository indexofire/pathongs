# CRAM 数据格式

---

CRAM文件格式是一种更密集的BAM文件形式，具有节省大量磁盘空间的优点。 虽然BAM文件包含文件中的所有序列数据，但CRAM文件通过利用额外的外部`Ref` - 参考序列文件而变小。 压缩和解压缩读取信息都需要此文件。

CRAM数据主要用于人等高等生物数据储存，微生物基因组比较小，一般不需要储存为CRAM格式。查看[此处](https://genome.ucsc.edu/goldenPath/help/cram.html)获取有关CRAM格式的更多信息。
