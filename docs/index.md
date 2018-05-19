# >写在开头

高通量测序技术已经发展了很多年，广受欢迎的小型测序仪[MiSeq][]和[PGM][]已经经历住市场的考验。同时新的技术和仪器的出现仍在不停推动测序技术领域的应用，如MinIon测序仪的现场实时测序；如中疾控对华大自主测序仪的应用等等。

随着Pulsenet China的逐渐完善，利用NGS技术对病原进行研究的工作迟早会成为一种CDCers的日常活动。这个工作包含了核酸提取、文库制备、测序、数据分析。前3项对于疾控中心实验室人员来说并非难事，但是数据分析对于医学或者生物学背景的CDCers来说往往成了最大的痛点。对于科研工作，可以与专业的生信公司人员合作。但是对于日常工作，目前总体上还是缺乏完善的工作流。即使是部署上GUI程序如[CLC Workbench](https://www.qiagenbioinformatics.com/products/clc-main-workbench/)，或者是类似[Galaxy](https://usegalaxy.org/)这样的web平台的应用，如果缺乏对数据分析最基本的概念，那也只能是空有一身装备却无法发挥。

杭州市疾控中心着眼于NGS的发展，早在2011年就开始进行有计划的学习和进修。经过长期的调研和准备，终于在2016年采购了[illumina][] MiSeq测序仪。虽然进行了长期的知识储备和积累，但是我们的工作人员也同样存在着对数据分析的难点。对于基层疾控中心实验室，受限于人员编制等问题，这些困难就只能由我们自己解决。所以也就有了这部教学内容。

笔者毕竟不擅长教学，因此内容和结构比较混乱，只能取名为心得笔记。如果对于刚进入这个领域的同行如能有一点点帮助，也算是了却心愿了。

## 版本更新

v0.0.3: 2018.05.15

 * 重新开始更新

v0.0.2: 2017.07.20

 * edirect 部分基本完成，进入矫错和补充状态
 * assembly 部分完成度50%

v0.0.1: 2017.04.01

 * 对原有内容做了适当调整
 * 重新调整了章节

[Linux]: http://www.linux.com/ "Linux"
[illumina]: http://www.illumina.com/ "Illumina"
[MiSeq]: http://www.illumina.com/search.ilmn?search=MiSeq&Pg=1&ilmn_search_btn.x=1 "MiSeq"
[PGM]: https://www.thermofisher.com/cn/zh/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-run-sequence/ion-pgm-system-for-next-generation-sequencing.html "PGM"
[Github]: https://www.github.com/ "Github"
