#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import re
import pysam
from pyfaidx import Fasta
import multiprocessing

from ddtools.utilities import *


def fecth_modif_site_wrapper(args):
	return fecth_modif_site(*args)

def fecth_modif_site(contig, bases=None):

	bf = pysam.AlignmentFile(bamfile, 'rb')
	bamRecords = bf.fetch(contig)

	for read in bamRecords:
		
		chrom = read.reference_name
		score = read.mapping_quality
		# filter reads
		
		if read.is_unmapped:
			continue
		if read.is_secondary:
			continue
		if score < mapq:
			continue

		# for PE:
		#		additional filter steps
		if read.is_paired:

			if read.is_read2:
				continue
			if read.mate_is_unmapped:
				continue
			if not read.is_proper_pair:
				continue
		
		
		# for negtive strand
		# predicted DNA damage site is locate in read end + 1
		if read.is_reverse:
			# removing reads contains adapter
			if re.match(".*\d+S$", read.cigarstring):
				continue

			# obtain the damage site
			start = read.reference_end
			end = read.reference_end + 1
			strand = "-"

			sequence = fetch_sequence(ref_fa, chrom, start, end)
			if sequence == None:
				continue
			seq_up = sequence.reverse.complement.seq.upper()
			#score = read.mapping_quality
		

		# for positive strand
		# predicted DNA damage site is locate in read start - 1
		else:
			if re.match("^\d+S.*", read.cigarstring):
				continue

			# obtain the damage site
			start = read.reference_start - 1
			end = read.reference_start
			strand = "+"

			sequence = fetch_sequence(ref_fa, chrom, start, end)
			if sequence == None:
				continue
			seq_up = sequence.seq.upper()
			#score = read.mapping_quality

		if len(seq_up) != 1:
			continue

		if seq_up not in ['A','C','G','T']:
			continue
		
		# output raw bed
		fo_raw_bed.write("\t".join(map(str, [chrom, start, end, seq_up, score, strand])) + '\n')

		'''
		for i in range(11):
			basetype = 'base_' + seq_up[i]
			raw_base_count.loc[i][basetype] += 1

		# if set filter
		if bases != None:

			median_base = seq_up[5]
			
			if median_base in bases:
				continue
			
			else:
				fo_fl_bed.write("\t".join(map(str, [chrom, start, end, seq_up, score, strand])) + '\n')
				
				for i in range(11):
					basetype = 'base_' + seq_up[i]
					fl_base_count.loc[i][basetype] += 1
		'''

	bf.close()

	return


	
def main(args=None):

	# args:
	# 	bam, fasta, out, select_base, mapq, plot,threads --deprecated
	#	bam, fasta, out, mapq, threads

	#bamfile = pysam.AlignmentFile(transferPath(args.bam),'rb')
	creat_raw_bed = FileIO(str(args.out), None, 'file') # a full path to file (create a new file)

	threadsUsing = args.threads


	# define global variables
	global bamfile, mapq, ref_fa, fo_raw_bed
	bamfile = transferPath(args.bam)
	bamfile_name = os.path.basename(bamfile)
	ref_fa = Fasta(transferPath(args.fasta))
	mapq = args.mapq

	fo_raw_bed = open(creat_raw_bed, 'w')

	'''
	deprecated

	if filter_flag:
		to_keep = set([x.upper() for x in args.select_base])
		
		# check if base type is valid
		diff = list(to_keep - set(['A','C','G','T']))
		if len(diff) > 0:
			info('{} were not valid bases.'.format(
				','.join(map(str, diff))
			), 'error')
			sys.exit()
		
		del diff

		to_keep = list(to_keep) # ['G']or ['A','G']
		to_keep_suffix = ''.join(map(str, to_keep)) # 'G' or 'AG'

		to_filter = [x for x in ['A','C','G','T'] if x not in to_keep]

		# create second new file to save filtered results
		fl_tmp = creat_raw_bed.split('.')[:-1]
		fl_tmp.insert(-1, to_keep_suffix)
		################################
		################################
		# need fixed
		#  test.bed.temp  test.bed.G.temp.temp
		################################
		################################
		creat_fl_bed = FileIO('.'.join(map(str, fl_tmp)) +'.temp', None, 'file')# full path
		global fo_fl_bed
		fo_fl_bed = open(creat_fl_bed, 'w')
	

	# check species
	species = 'human'
	fasta_chroms = sorted([x for x in ref_fa.keys() if '_' not in x])
	query = set(['chr20','chr21','chr22'])
	if len(query - set(fasta_chroms)) == 3:
		species = 'mouse'

	contigs = ['chr' + str(i) for i in range(1,20)] + ['chrX']
	if species == 'human':
		contigs.extend(['chr20','chr21','chr22'])

	contigs.sort()


	info('Presuming the species is {} according to the Fasta file you provided.'.format(species), 'out')
	
	############
	# test
	#print('\nraw_bed_path')
	#print(creat_raw_bed)
	#print('\nfl_bed_path')
	#print(creat_fl_bed)
	#print('\nchromosomes')
	#print(contigs)
	############
	'''
	contigs = ref_fa.keys()
	info('Start allocate jobs.', 'out')
	TASKS = []
	# split bam into chunks by chromosome name
	for contig in contigs:
		TASKS.append(tuple([contig]))
	
	# allocate to multiple cores
	info('Start running jobs.', 'out')
	pool = multiprocessing.Pool(threadsUsing)
	pool.map_async(fecth_modif_site_wrapper, TASKS)
	pool.close()
	pool.join()

	'''
	Error:
	TypeError: self.b,self.htsfile,self.index,self.iter cannot be converted to a Python object for pickling

		已经open的文件，无法通过参数传递， 因为无法 “序列化”， 检查一个变量是否可以序列化的方法是：

				import pickle
				var = a
				pickle.dumps(var)

	'''


	# close bed files and sort bed using UCSC tools
	fo_raw_bed.close()
	info('{} covert to Bed Done.'.format(bamfile_name), 'out')


	'''
	deprecated:
		this function only represents convert format & filter

	sort_cmd = "bedSort {} {} && rm {}".format(creat_raw_bed, creat_raw_bed[:-5], creat_raw_bed)
	p = os.system(sort_cmd)
	if filter_flag:
		fo_fl_bed.close()
		sort_cmd = "bedSort {} {} && rm {}".format(creat_fl_bed, creat_fl_bed[:-5], creat_fl_bed)
		p = os.system(sort_cmd)

	# output base counts
	# df.to_csv(fo,sep='\t',float_format=None,header =True,index=True, quoting=None)
	mdf = []
	creat_fo_bases = FileIO(
		'.'.join(map(str, creat_raw_bed[:-5].split('.')[:-1]+['bases'])), None, 'file')
	fo_bases = open(creat_fo_bases, 'w')
	fo_bases.write('BeforeFilter:\n')
	raw_base_count.to_csv(fo_bases, sep='\t', float_format=None, header =True, index=True, quoting=None)
	mdf.append(raw_base_count)

	if filter_flag:
		fo_bases.write('\nAfterFilter:\n')
		fl_base_count.to_csv(fo_bases, sep='\t', float_format=None, header =True, index=True, quoting=None)
		mdf.append(fl_base_count)
	
	fo_bases.close()

	# close genome files
	ref_fa.close()

	info('output done. start plot', 'out')

	# plot base composition
	if args.plot:
		figname = '.'.join(map(str, creat_raw_bed[:-5].split('.')[:-1]+['bases', 'pdf']))
		plotBase(mdf, figname)
	'''

if __name__ == '__main__':
	main()
