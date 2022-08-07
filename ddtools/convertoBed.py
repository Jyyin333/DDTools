#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import re
import time
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
	info('Loading data...', 'out')
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

	# close bed files and sort bed using UCSC tools
	fo_raw_bed.close()
	
	info('Start sort bed and build index.', 'out')
	#####################
	# 4.20 adding
	# sort and idnex
	#####################
	_temp_bed = creat_raw_bed + '.' + get_random_string(8)
	bedsort_cmd = 'bedSort {} {}'.format(creat_raw_bed, _temp_bed)

	# sort Bed
	os.system(bedsort_cmd)
	latest_size = 1
	zipfile_size = -1
	while latest_size > zipfile_size:
		zipfile_size = latest_size
		latest_size = os.path.getsize(_temp_bed)
		time.sleep(2)

	time.sleep(2)
	del latest_size, zipfile_size


	# rename and bgzip
	os.remove(creat_raw_bed)
	os.rename(_temp_bed, creat_raw_bed)

	bgzip_cmd = 'bgzip {}'.format(creat_raw_bed)
	bgziped_file = creat_raw_bed + '.gz'
	tabix_cmd = 'tabix -p bed {}'.format(bgziped_file)

	os.system(bgzip_cmd)
	latest_size = 1
	zipfile_size = -1
	while latest_size > zipfile_size:
		zipfile_size = latest_size
		latest_size = os.path.getsize(bgziped_file)
		time.sleep(2)

	time.sleep(2)
	os.system(tabix_cmd)
	time.sleep(10)

	info('{} covert to Bed Done.'.format(bamfile_name), 'out')


if __name__ == '__main__':
	main()
