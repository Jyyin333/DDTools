# DDTools使用说明

DDTools是一款基于Python编写的、用于结合末端标记法策略的二代测序相关数据的分析工具  
包含项目构建、数据质控、数据分析和图形展示等多个功能模块，可以帮助用户实现对数据的快速处理

![DDTools modules](.\figures\ddtools_modules.png)

## DDTools的下载与安装
目前，用户可以通过下载本页面中的源码，通过```python setup.py install```安装使用  
未来，也计划将DDTools上传至PyPi，届时可通过pip的方式直接下载

## DDTools presentation
接下来以一套CLAPS-seq的数据为例，演示DDTools的使用方法及流程，原始测序数据来源于[GSE181312](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181312)  
由于文件大小限制问题，本页面无法将原始及中间文件进行分享，请按照个人情况自行下载

| FILE | DESCRIPTION |
| :---: | :---: |
| GSM5494639 | HeLa-S3 INPUT rep1 |
| GSM5494647 | HeLa-S3 VIVO rep1 |

## 软件准备
用户需要提前安装一些必备软件：

+ [bedSort](http://hgdownload.cse.ucsc.edu/admin/exe/)
+ htslib

此外，DDTools中包含的**makeProject**命令，可以按照用户指定的参数，帮助建立[snakemake](https://snakemake.readthedocs.io/en/stable/)流程脚本完成上游分析  
其中会涉及到多个对其它工具的使用，因此，使用该命令需要做一些额外的准备，即提前在shell环境中安装更多的软件  
所需软件如下：

+ bwa / bowtie / bowtie2
+ picard
+ fastqc
+ samtools
+ cutadapt
+ snakemake

> 本节演示主要从**makeProject**命令开始，因此推荐以FASTQ作为输入文件
> 如果用户是从SAM/BAM格式的比对文件开始，此步骤可以跳过，并使用**convertoBed**命令进行处理
