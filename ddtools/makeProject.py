#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import sys, os
import argparse
import json
import gzip
import time
from collections import defaultdict

from ddtools.utilities import *


def searh_fastq_files(target):

	'''
	return:
		a list a absolute paths of fastq file
		['/home/usr/fastq/xxx.fastq.gz']
	'''

	#res = {'samples':{}}
	res = [] # store path to fastq

	try:
		assert os.path.isdir(target)
	except AssertionError:
		info('The {} is not a directory.'.format(target), 'error')
		sys.exit()

	cwd_files = os.listdir(target)
	non_dirs = [x for x in cwd_files if os.path.isfile(os.path.join(target,x))]
	sub_dirs = [x for x in cwd_files if os.path.isdir(os.path.join(target,x))]

	if not len(sub_dirs) == 0:
		# Probably because each sample has its own folder

		for sub_dir in sub_dirs:

			dirname = os.path.join(target,sub_dir)
			prob_fq_files = [os.path.join(dirname, x) for x in os.listdir(dirname) if x.endswith('.fq.gz') or x.endswith('.fastq.gz') or x.endswith('.fastq') or x.endswith('.fq')]

			res.extend(prob_fq_files)


	prob_fq_files = [os.path.join(target, x) for x in non_dirs if x.endswith('.fq.gz') or x.endswith('.fastq.gz') or x.endswith('.fastq') or x.endswith('.fq')]
	res.extend(prob_fq_files)

	# all search done.
	return res



def makeSampleList(dirname, namelist):

	'''
	# dirname:
	# 		/usr/home/project/fastq_tmp
	#
	# an example:
	# 		R21022539-snOGseq8-Oggi6_combined_R1.fastq.gz

	returns a dict:
		collect = {
			'samples':{
				'R21022539-snOGseq8-Oggi6_combined':{
					'id':xx,
					'pl':xx,
					'pu':xx
				}
			},

			'fastq':{
				'R21022539-snOGseq8-Oggi6_combined':{
					'Read1':abspath2fq,
					'Read2':abspath2fq,
			}
		}
	'''

	collect = {'samples':{}, 'fastq':{}}
	samples = []
	fqfiles = os.listdir(dirname)
	
	if len(fqfiles):
		#info('Fastq files were detected under {}, stop searching for secondary folder.'.format(target), 'warn')
		suffixs = set()
		for sm in fqfiles:
			line = sm.split('.')
			sample = line[0]
			sf = '.'.join(map(str, line[1:]))
			samples.append(sample)
			suffixs.add(sf)

		if not len(suffixs) == 1:
			info('Fastq files are not consistent in suffixs.', 'error')
			sys.exit()
		
		# for single-end sequencing
		test = [samples[0].rfind('_R1'), samples[0].rfind('_R2')]
		if test.count(-1) == 2:
			info('Not support single-end sequencing data.', 'error')
			sys.exit()
		
		sample_res = defaultdict(list)
		for line in samples:
			if line.rfind('_R1') != -1:
				sample_base = line[:line.rfind('_R1')]
				reads_group = line[line.rfind('_R1'):]
				sample_res[sample_base].append(reads_group)
				sample_res[sample_base].append(reads_group.replace('_R1','_R2',1))
	
	# collect sample info
	suffix = '.' + str(list(suffixs)[0])
	for sm,values in sample_res.items():
		# sm = R21022539-snOGseq8-Oggi6
		# values = ['_R1','_R2']
		collect['fastq'][sm] = {}
		collect['fastq'][sm]['Read1'] = os.readlink(os.path.join(dirname, sm + values[0] + suffix))
		collect['fastq'][sm]['Read2'] = os.readlink(os.path.join(dirname, sm + values[1] + suffix))

		# get pu using R1 fastq file
		if not suffix.endswith('.gz'):
			with open(collect['fastq'][sm]['Read1']) as fi:
				index = fi.readline().rstrip().split(':')[-1]
		else:
			with gzip.open(collect['fastq'][sm]['Read1']) as fi:
				index = fi.readline().decode().rstrip().split(':')[-1]

	########################
	# 4.20
	# to be modified
	# if fastq files are downloaded from GEO
	# some keys within this dict may not applicable
	########################

		collect['samples'][sm] = {}
		collect['samples'][sm]['id'] = sm
		collect['samples'][sm]['lb'] = namelist.get('library')
		collect['samples'][sm]['pl'] = namelist.get('platform')
		collect['samples'][sm]['pu'] = index
		collect['samples'][sm]['sample'] = sm
	
	return collect



def main(args=None):

	info('Start makeProject.')
	cwd = os.getcwd()

	snake_threads = args.snake_threads
	snakeRun = args.run

	namelist = {}
	namelist['project_name'] = args.project_name
	if namelist['project_name'] == None:
		namelist['project_name'] = os.path.basename(cwd)
		
	namelist['platform'] = args.platform
	namelist['library'] = args.library_name
	namelist['aligner'] = args.aligner

	genomelist = {}
	genomelist['fasta'] = transferPath(args.fasta)
	genomelist['index'] = transferPath(args.aligner_index)

	# check if genome files are exist
	if not os.path.exists(genomelist['fasta']):
		info('{} not exists.'.format(genomelist['fasta']), 'error')
		sys.exit()

	if not checkIndex(genomelist['index']):
		info('{} is not a valid {} index path.'.format(genomelist['index'], namelist['aligner']), 'error')
		sys.exit()

	# adding picard path
	# genomelist['picard'] = transferPath(args.picard)
	while True:

		# check if picard exist in environment variable
		if checkUnixTools('picard'):
			picard_path = transferPath(os.popen('which picard').read().strip())
			path_pre = '/'.join(map(str,picard_path.split('/')[:-2]))
			_cmd = 'find {} -name picard.jar'.format(path_pre)
			picardJar = os.popen(_cmd).read().strip()
			time.sleep(10)
			if len(picardJar) == 0:
				picardJar = None
			genomelist['picard'] = picardJar
			break

		else:
			genomelist['picard'] = None
			break

	if genomelist['picard'] == None:
		if args.picard:
			genomelist['picard'] = transferPath(args.picard)
		else:
			info('No picard in your environment path. Please provide the path to picard.jar', 'error')
			sys.exit()


	# 1/12 adding
	genomelist['mapq'] = args.mapq
	'''
	selectbases = []
	selectbases_line = ''
	for i in args.select_base:
		i = str(i).upper()
		selectbases.append(i)
		selectbases_line += i
		selectbases_line += ' '
	
	genomelist['select_base'] = selectbases_line.rstrip()
	'''

	# check whether the current path is correct
	cwd_files = os.listdir()
	if 'fq' not in cwd_files and 'fastq' not in cwd_files:
		info('Current work dirctory do not contains ./fastq or ./fq dir. Please check the path.\n\tNote: Use lowercase letters', 'error')
		sys.exit()

	try:
		assert cwd_files.count('fq') + cwd_files.count('fastq') == 1
	except AssertionError:
		info('Seemingly you have two directories for fastq files, please put them all in one folder', 'error')
		sys.exit()


	# search for fastq files
	fq_dir = 'fq' if 'fq' in cwd_files else 'fastq'
	fq_dir = os.path.join(cwd, fq_dir)
	pfq = searh_fastq_files(fq_dir)
	if len(pfq) == 0:
		info('Do not find fastq files in {}. Please check your directory.'.format(fq_dir), 'error')
		sys.exit()
		
	# create temp dir to link fastq files
	new_fq_dir = FileIO(cwd, 'temp_fastq', 'dir')

	for fqpath in pfq:
		name = os.path.basename(fqpath)
		os.symlink(fqpath, os.path.join(new_fq_dir,name))

	# softlink done.

	# make sample list and json file
	res = makeSampleList(new_fq_dir, namelist)
	sampleList = sorted(res.get('samples').keys())

	# export file
	initial_dir = FileIO(cwd, 'initial', 'dir')

	sname = namelist['project_name'] + '.sampleList'
	with open(os.path.join(initial_dir, sname), 'w') as smout:
		for i in sampleList:
			smout.write(i+'\n')
	
	pname = namelist['project_name'] + '_project.json'
	with open(os.path.join(initial_dir, pname), 'w') as jsonout:
		json.dump(res, jsonout, indent=4)

	# remove fastq temp dir
	os.chdir(new_fq_dir)
	for i in os.listdir():
		os.remove(i)
	os.chdir(cwd)
	os.rmdir(new_fq_dir)



	# Construct snakefile
	info('Construct snakefile.')
	snakefile = FileIO(cwd, namelist['project_name']+'.snakefile', 'file')
	mysnake = SnakeBlock(namelist['project_name'])
	mysnake.constructVar(genomelist) # add genome path to snakefile
	
	## add steps
	### fastqc
	fastqc_list = [
	'input',
	'r1=lambda wildcards:config["fastq"][wildcards.sample]["Read1"]',
	'r2=lambda wildcards:config["fastq"][wildcards.sample]["Read2"]',
	'output',
	'touch("dataQC/Sample_{sample}/{sample}.fastQC.Done")',
	'params',
	'outDir="dataQC/Sample_{sample}"',
	'threads',
	'2',
	'shell',
	'"fastqc --noextract --format=fastq -q --threads={threads} --outdir {params.outDir} {input.r1} {input.r2}"'
	]
	mysnake.constructBlock('fastQC', fastqc_list, collect=True)

	### preprocessing
	if args.cutadapt_params:
		cutadapt_list = [
		'input'
		]
		pass

	### mapping to genome
	if namelist['aligner'] == 'bwa':
		map_list = [
		'input',
		'r1=lambda wildcards:config["fastq"][wildcards.sample]["Read1"]',
		'r2=lambda wildcards:config["fastq"][wildcards.sample]["Read2"]',
		'output',
		'"mapped/{sample}.bam"',
		'params',
		'rg="@RG\\\\tID:{sample}\\\\tSM:{sample}"',
		'lb=lambda wildcards:config["samples"][wildcards.sample]["lb"]',
		'pl=lambda wildcards:config["samples"][wildcards.sample]["pl"]',
		'pu=lambda wildcards:config["samples"][wildcards.sample]["pu"]',
		'threads',
		'5',
		'shell',
		'''"bwa mem -t {threads} -R '{params.rg}\\\\tLB:{params.lb}\\\\tPL:{params.pl}\\\\tPU:{params.pu}' {ALIGNER_INDEX} {input.r1} {input.r2} | samtools view -Sb - > {output}"'''
		]

	if namelist['aligner'] == 'bowtie':
		pass

	if namelist['aligner'] == 'bowtie2':
		map_list =[
		'input',
		'r1=lambda wildcards:config["fastq"][wildcards.sample]["Read1"]',
		'r2=lambda wildcards:config["fastq"][wildcards.sample]["Read2"]',
		'output',
		'"mapped/{sample}.bam"',
		'threads',
		'5',
		'shell',
		'''"bowtie2 --no-mixed --no-discordant -p {threads} -x {ALIGNER_INDEX} -1 {input.r1} -2 {input.r2} | samtools view -Sb - > {output}"'''
		]

	mysnake.constructBlock('mapping', map_list)

	### sort bamfile
	sortbam_lsit = [
	'input',
	'rules.mapping.output',
	'output',
	'bam=temp("mapped/{sample}_sorted.bam")',
	'bai=temp("mapped/{sample}_sorted.bai")',
	'threads',
	'5',
	'shell',
	'''"java -Xmx40g -jar {picard} SortSam -I {input} -O {output.bam} --SORT_ORDER coordinate --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"'''
	]
	mysnake.constructBlock('SortBam', sortbam_lsit)

	### remove duplicates
	dedup_list = [
	'input',
	'bam="mapped/{sample}_sorted.bam"',
	'bai="mapped/{sample}_sorted.bai"',
	'output',
	'bam="mapped/{sample}/{sample}_sorted_dedup.bam"',
	'bai="mapped/{sample}/{sample}_sorted_dedup.bai"',
	'metrics="mapped/{sample}/{sample}_sorted_dedup.metrics"',
	'threads',
	'5',
	'shell',
	'''"java -Xmx40g -jar {picard} MarkDuplicates -I {input.bam} -O {output.bam} --METRICS_FILE {output.metrics} --REMOVE_DUPLICATES true --ASSUME_SORTED true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"'''
	]
	mysnake.constructBlock('Dedup', dedup_list)

	### collectMultipleMetrics
	cmm_list = [
	'input',
	'bam="mapped/{sample}_sorted.bam"',
	'bai="mapped/{sample}_sorted.bai"',
	'output',
	'"dataQC/Sample_{sample}/{sample}.alignment_summary_metrics"',
	'params',
	'prefix="dataQC/Sample_{sample}/{sample}"',
	'threads',
	'5',
	'shell',
	'''"java -Xmx20g -jar {picard} CollectMultipleMetrics -I {input.bam} -O {params.prefix} -R {GENOME_FA}"'''
	]
	mysnake.constructBlock('collectMultipleMetrics', cmm_list, collect=True)

	### fetch predicted modification site
	flist = [
	'input',
	'bam="mapped/{sample}/{sample}_sorted_dedup.bam"',
	'bai="mapped/{sample}/{sample}_sorted_dedup.bai"',
	'output',
	'fbed="beds/{sample}/{sample}.clean.bed"',
	'threads',
	'1',
	'shell',
	'''"dmtools convertoBed -b {input.bam} -g {GENOME_FA} -o {output} -q {MAPQ} -t {threads}"'''
	]
	mysnake.constructBlock('fetch_modification_sites', flist)


	### bedtobam
	tobamlist = [
	'input',
	'rules.fetch_modification_sites.output.fbed',
	'output',
	'"mapped/{sample}/{sample}.clean.bam"',
	'shell',
	'''"bedtools bedtobam -g {GENOME_FA} -mapq {MAPQ} -i {input} > {output} && sleep 1s && samtools index {output}"'''
	]
	mysnake.constructBlock('bedtobam', tobamlist)

	if checkUnixTools('deeptools'):
		### bamtobigwig
		tobwlist = [
		'input',
		'rules.bedtobam.output',
		'output',
		'"bigwigs/{sample}.rpkm.bw"',
		'threads',
		'5',
		'shell',
		'''"bamCoverage -p {threads} --normalizeUsing RPKM -b {input} -of bigwig -o {output} -bs 10"'''
		]
		mysnake.constructBlock('bamtobigwig', tobwlist, collect=True)

	# make done & write out
	mysnake.constructHeadAndTail()
	mysnake.writeFile(snakefile)

	info('Generate DAG.')
	# generate DAG
	dryrun = 'snakemake -s {} --dag |dot -Tsvg > {}'.format(
		snakefile, snakefile+'.svg')
	os.system(dryrun)
	time.sleep(10)

	# run snakefile
	if snakeRun:
		runsnakefile_cmd = 'nohup snakemake -s {} --latency-wait 5 -j {} > {} 2>&1 &'.format(snakefile, snake_threads, snakefile+'.running.log')
		os.system(runsnakefile_cmd)


	# Done
	info('makeProject Done.')

if __name__ == '__main__':
	main()

