# DDTools使用说明

DDTools是一款基于Python编写的、用于结合末端标记法策略的二代测序相关数据的分析工具  
包含项目构建、数据质控、数据分析和图形展示等多个功能模块，可以帮助用户实现对数据的快速处理

![DDTools modules](./figures/ddtools_modules.png)

## DDTools的下载与安装
目前，用户可以通过下载本页面中的源码，通过```python setup.py install```安装使用  
未来，也计划将DDTools上传至PyPi，届时可通过pip的方式直接下载

## DDTools presentation
接下来以一套CLAPS-seq的数据为例，演示DDTools的使用方法及流程，ClAPS-seq是一种用于在全基因组范围内检测8-oxodG氧化损伤的方法，原始测序数据来源于[GSE181312](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181312)  
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

> 推荐以FASTQ作为输入文件，完整使用DDTools进行数据的全套分析
> 如果用户是从SAM/BAM格式的比对文件开始，此步骤可以跳过，并使用**convertoBed**命令进行处理

## DDTools Workflow
现在，我们有两个样本，将其放在fastq/目录下，并使用**makeProject**命令进行前期处理
```
# pwd: /fastq/
.
├── claps-seq-vivo-rep1
│	├── claps-seq-vivo-rep1_R1.fastq.gz
│	└── claps-seq-vivo-rep1_R2.fastq.gz
│
└── claps-seq-input
	├── claps-seq-input_R1.fastq.gz
	└── claps-seq-input_R2.fastq.gz
```
```bash
ddtools makeProject --project_name proj-test \
--library_name claps-seq \
--fasta $path_to_genome_fasta \
--aligner bwa --aligner_index $path_to_aligner_index \
--mapq 25 \
--snake_threads 6 --run
```
通过上述命令，即可完成数据的上游分析，即数据的质控、比对、过滤等  
用户可以将参数保存到一个文件中，便于以后针对相同类型数据使用，命令如下
```bash
# 保存参数到params.list
cat << EOF > params.list 
--project_name
proj-test
--library_name
claps-seq
--fasta
path_to_genome_fasta
--aligner
bwa
--aligner_index
path_to_aligner_index
--mapq
25
--snake_threads
6
--run
EOF
# 运行makeProject，效果同上
ddtools makeProject @params.txt
```

在参数中，通过--project_name给此次项目命名，后续产生的中间文件会以该名字作为前缀；--aligner告诉DDTools使用哪一种工具进行比对，目前支持bwa、bowtie和bowtie2；  
--mapq设置最小比对质量用于对reads进行筛选；--snake_threads设置cpu核心的使用数目；--run则是告诉DDTools即刻运行创建的snakefile文件
> 更为详细的参数介绍，可以通过--help查看帮助文档

运行完成后，会得到经过初步过滤的BED格式的文件，该文件记录了理论上氧化损伤在基因组上的出现位置

如果用户从SAM/BAM文件开始，则可以通过**convertoBed**对数据进行过滤和格式转换操作，同样可以得到适用于DDTools下游的BED文件
```bash
ddtools convertoBed --bam input.bam --fasta genome.fa --mapq 25
```

得到BED文件后，便可以通过DDTools的不同模块进行数据的探索  
在DDTools中，目前包含了四个下游分析模块，分别是seqContext、regAnnotation、computeGbias和computeMTX以及配套的绘图命令

1. seqContext    用于分析理论损伤位点处及附近碱基序列组成
2. regAnnotation    用于损伤位点的注释
3. computeGbias    在更大的局部范围内查看损伤位点分布是否有序列偏好性
4. computeMTX    用于查看在某一类型基因组区域（如transcription start site）的损伤分布情况
