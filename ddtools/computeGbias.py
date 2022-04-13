#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import pysam, gzip
import numpy as np
import multiprocessing
from pyfaidx import Fasta
import py2bit

from ddtools.utilities import *


def _fetchSizefromtbit(tbit):
	'''
	returns:
		(INT)
		Effective Genome Size retrived from genome 2bit file.

		(List of list)
		Genome Chunks
	'''
	chroms = tbit.chroms()
	common_chroms = [x for x in chroms.keys() if not x.count('_')]
	res = 0
	GenomeChunks = []
	for i in common_chroms:
		res += chroms.get(i)

	for chrom in common_chroms:
		size = chroms.get(chrom)
		start_pos = 0
		for start in np.arange(start_pos, size, 1e5):
			end_pos = min(start+1e5, size)

			GenomeChunks.append([chrom, start, end_pos])

	return res, GenomeChunks



def computeG_wrapper(args):
	return computeG_worker(*args)

def computeG_worker(chrom, start, end, stepsize, regionsize, treatGroups, tbit_file, controlFile=None):
	
	treatBedfiles = [pysam.TabixFile(bed, mode='r', parser=pysam.asBed()) for bed in treatGroups]
	tbit = py2bit.open(tbit_file)

	sub_signal_per_gc = []
	positions = np.arange(start, end, stepsize)

	if controlFile:
		ctbx = pysam.TabixFile(controlFile, mode='r', parser=pysam.asBed())

	for idx in range(len(positions)):
		frag_start = int(positions[idx])
		## 2022/3/25 adding
		## in case of error
		## RuntimeError: you must supply a chromsome
		frag_end = frag_start + regionsize
		if end < frag_end:
			break

		# fetch G/C contents
		base_G = get_base_content(chrom, frag_start, frag_end, tbit, fraction=False, base='G')
		base_C = get_base_content(chrom, frag_start, frag_end, tbit, fraction=False, base='C')

		# fetch read counts
		treatCounts = []
		for idx, tbx in enumerate(treatBedfiles):
			tbx_p, tbx_m = bedCountsinRegion(tbx, chrom, frag_start, frag_end)
			treatCounts.extend([tbx_p, tbx_m])

		# fetch control counts if exist
		controlCounts = []
		if controlFile:
			ctbx_p, ctbx_m = bedCountsinRegion(ctbx, chrom, frag_start, frag_end)
			controlCounts.extend([ctbx_p, ctbx_m])

		treatCounts.extend(controlCounts)
		sub_signal_per_gc.append(tuple([*treatCounts, base_G, base_C]))


	tbit.close()

	for t in treatBedfiles:
		t.close()

	if controlFile:
		ctbx.close()

	return sub_signal_per_gc



def main(args=None):

	# loading params
	samplesize = args.samplesize
	regionsize = args.regionsize
	cores = args.cores
	tbit_file = transferPath(args.genome2bit)
	outfile = transferPath(args.output)

	## check bedfile is valid
	bedfiles = [transferPath(x) for x in args.bed]
	nSamples = len(bedfiles)
	controlBed = args.control

	for bed in bedfiles:
		try:
			tbx = pysam.TabixFile(bed, mode='r', parser=pysam.asBed())
			tbx.close()
		except OSError as oe:
			info()
			sys.exit()

	if controlBed:
		controlBed = transferPath(controlBed)
		try:
			tbx = pysam.TabixFile(controlBed, mode='r', parser=pysam.asBed())
			tbx.close()
		except OSError as oe:
			info()
			sys.exit()


	# split genome into chunks
	tbit = py2bit.open(tbit_file)
	effectiveSize, GenomeChunks = _fetchSizefromtbit(tbit)
	stepsize = int(effectiveSize / samplesize)
	tbit.close()


	# generate tmp bed file
	ref_fa = Fasta(transferPath(args.genomefasta))
	keep_bases = args.keep_bases
	# tmp_bedfiles: a list of new bedfiles name with suffix of random string
	tmp_bedfiles = [str(bed) + '.' + get_random_string(8) for bed in bedfiles]
	# tmp_bedfiles_list: a list of lists of new bedfiles name with suffix of .gz/.tbi
	tmp_bedfiles_list = [filterBed(bedfiles[i], ref_fa, tmp_bedfiles[i], keep_bases) for i in range(nSamples)]
	ref_fa.close()

	inputBeds = [x[0] for x in tmp_bedfiles_list]
	# running main procedure
	JOBS = []
	for reg in GenomeChunks:
		arglist = [reg[0], int(reg[1]), int(reg[2]), stepsize, regionsize, inputBeds, tbit_file, controlBed]
		JOBS.append(tuple(arglist))

	pool = multiprocessing.Pool(cores)
	res = pool.map_async(computeG_wrapper, JOBS).get()
	pool.close()
	pool.join()


	# create headerline
	headerLine = {
		'Bedfiles':bedfiles,
		'ControlFile':controlBed,
		'SampleSize':samplesize,
		'RegionSize':regionsize,
		'GenomeFasta':str(transferPath(args.genomefasta)),
		'Genome2bit':str(transferPath(args.genome2bit)),
		'keep_bases':keep_bases,
		'OutputFile':str(outfile)
		}

	writeHeader(outfile, headerLine)


	with open(outfile, 'a') as fo:
		
		for signal in res:
			
			for line in signal:
				line = '\t'.join(map(str, line)) + '\n'
				fo.write(line)

	# remove tmp bedfile
	for f in tmp_bedfiles_list:
		os.remove(f[0])
		os.remove(f[1])

if __name__ == '__main__':
	main()

