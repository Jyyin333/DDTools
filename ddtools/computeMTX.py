#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import math, gzip
import pysam
from pyfaidx import Fasta
import numpy as np
import pandas as pd
import multiprocessing

from ddtools.utilities import *


def fetchEffectiveGenomeSize(genomebuild):
	
	effectivegenomesize = {
	'grch37':2864785220,
	'hg19':2864785220,
	'grch38':2913022398,
	'hg38':2913022398,
	'mgscv37':2620345972,
	'mm9':2620345972,
	'grcm38':2620345972,
	'mm10':2620345972
	}

	genomebuild = genomebuild.lower()
	if genomebuild in effectivegenomesize.keys():
		return effectivegenomesize.get(genomebuild)

	return None


def refpoint_run_wrapper(args):
	return refpoint_run(*args)

def refpoint_run(bedfile, upstream, downstream, binsize, stepsize, totalCounts, method, r_l, g_s, sub_REGIONS):

	TS, NTS = [], []

	tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())

	for chrom, ref_site, ref_strand in sub_REGIONS:

		ts, nts = [], []

		if ref_strand == '+':

			region_up = ref_site - upstream
			region_down = ref_site + downstream

			#records = tbx.fetch(chrom, region_up, region_down)
			#records = list(records)

			sites = np.arange(region_up, region_down, stepsize, dtype=int)

			for idx in range(len(sites)):
				istart = sites[idx]
				iend = istart + binsize

				p, m = bedCountsinRegion(tbx, chrom, istart, iend, totalCounts, method, r_l, g_s) # plus/minus strand read counts

				ts.append(m)
				nts.append(p)

		elif ref_strand == '-':

			region_up = ref_site + upstream
			region_down = ref_site - downstream

			#records = tbx.fetch(chrom, region_up, region_down)
			#records = list(records)

			sites = np.arange(region_up, region_down, -stepsize, dtype=int)

			for idx in range(len(sites)):
				iend = sites[idx]
				istart = iend - binsize

				p, m = bedCountsinRegion(tbx, chrom, istart, iend, totalCounts, method, r_l, g_s) # plus/minus strand read counts

				ts.append(p)
				nts.append(m)

		else:
			continue

		TS.append(ts)
		NTS.append(nts)

	tbx.close()

	return TS, NTS



def refpoint_main(params_dict, REGION_LIST):

	params_dict = params_dict.copy()
	# Determine the type of Score file.
	if params_dict.get('genome2bit') != None:
		pass # return

	bedfile = params_dict.get('bedfile') # already converted into absolute path and could be opend by pysam.TabixFile
	tokeep = params_dict.get('keep_bases')
	fasta = params_dict.get('genomefasta')
	ref_fa = Fasta(fasta)

	#######
	# create tmp bed file to save subtracted reads
	# these tmp files would be deleted after this procedure completed.
	#######
	tmp_bedfile = str(bedfile) + '.' + get_random_string(8)
	tmpfile_list = filterBed(bedfile, ref_fa, tmp_bedfile, tokeep)
	ref_fa.close()
	#################### Done


	upstream = params_dict.get('upstream')
	downstream = params_dict.get('downstream')
	binsize = params_dict.get('binsize')
	stepsize = params_dict.get('stepsize')
	cores = params_dict.get('cores')
	totalCounts = params_dict.get('totalCounts')

	normlizeMethod = params_dict.get('norm_m')
	r_l = params_dict.get('readLength')
	g_s = params_dict.get('genomeSize')


	nbins = (int(abs(downstream + upstream)) // stepsize)
	params_dict['nBins'] = nbins

	df_ts = pd.DataFrame(0, index=[], columns=range(nbins), dtype=float)
	df_nts = pd.DataFrame(0, index=[], columns=range(nbins), dtype=float)

	# chunk REGION_LIST into small set of REGIONS respect to given cores
	nRegion = len(REGION_LIST)
	if nRegion <= cores:
		cores = nRegion
		nCut = 1
	else:
		nCut = math.ceil(nRegion / cores)


	JOBS = []
	i = 0
	while i < nRegion:
		
		sub_REGIONS = REGION_LIST[i:i+nCut]
		sub_job = tuple([tmpfile_list[0], upstream, downstream, binsize, stepsize, totalCounts, normlizeMethod, r_l, g_s, sub_REGIONS])
		JOBS.append(sub_job)
		i += nCut

	pool = multiprocessing.Pool(cores)
	res = pool.map_async(refpoint_run_wrapper, JOBS).get()
	pool.close()
	pool.join()

	# res : ([TS,NTS], [], [], ...)
	# TS/NTS: [[4,1,2,5,6], [5,3,2,6,7], [], ...]
	for TS, NTS in res:

		sub_TS = pd.DataFrame(np.asarray(TS, dtype=float).reshape(len(TS), nbins), columns=range(nbins))
		sub_NTS = pd.DataFrame(np.asarray(NTS, dtype=float).reshape(len(NTS), nbins), columns=range(nbins))

		# ensure the columns are same in length 
		#assert sub_TS.shape[1] == sub_NTS.shape[1] == nbins

		df_ts = pd.concat([df_ts,sub_TS], axis=0)
		df_nts = pd.concat([df_nts,sub_NTS], axis=0)

	# merge TS/NTS dataframe
	DF = pd.concat([df_ts, df_nts], axis=1)
	dim = DF.shape
	DF.index = range(dim[0])
	DF.columns = range(dim[1])

	# remove temp files
	#############################
	# Annotate for Testing
	for f in tmpfile_list:
		os.remove(f)
	#############################
	return DF, params_dict.copy()



def scale_run_wrapper(args):
	return scale_run(*args)

def scale_run(bedfile, up_bins, mid_bins, down_bins, totalCounts, method, r_l, g_s, sub_REGIONS):

	TS, NTS = [], []

	tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())

	for chrom, ref_strand, up_region, mid_region, down_region in sub_REGIONS:

		ts, nts = [], []

		if ref_strand == '+':

			# up-region
			up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
			up_pos = np.append(up_pos,up_region[1])
			for i in range(len(up_pos)-1):
				sub_start = up_pos[i]
				sub_end = up_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s) # plus/minus strand read counts
				ts.append(m)
				nts.append(p)

			# mid-region
			mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
			mid_pos = np.append(mid_pos,mid_region[1])
			for i in range(len(mid_pos)-1):
				sub_start = mid_pos[i]
				sub_end = mid_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s)
				ts.append(m)
				nts.append(p)

			# down-region
			down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
			down_pos = np.append(down_pos,down_region[1])
			for i in range(len(down_pos)-1):
				sub_start = down_pos[i]
				sub_end = down_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s)
				ts.append(m)
				nts.append(p)

		if ref_strand == '-':

			# down-region
			down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
			down_pos = np.append(down_pos,down_region[1])
			for i in range(len(down_pos)-1):
				sub_start = down_pos[i]
				sub_end = down_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s)
				ts.append(p)
				nts.append(m)

			# mid-region
			mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
			mid_pos = np.append(mid_pos,mid_region[1])
			for i in range(len(mid_pos)-1):
				sub_start = mid_pos[i]
				sub_end = mid_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s)
				ts.append(p)
				nts.append(m)

			# upstream-region
			up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
			up_pos = np.append(up_pos,up_region[1])
			for i in range(len(up_pos)-1):
				sub_start = up_pos[i]
				sub_end = up_pos[i+1]
				p, m = bedCountsinRegion(tbx, chrom, sub_start, sub_end, totalCounts, method, r_l, g_s)
				ts.append(p)
				nts.append(m)

			# reverse values
			ts = ts[::-1]
			nts = nts[::-1]

		else:
			continue

		TS.append(ts)
		NTS.append(nts)

	tbx.close()

	return TS, NTS



def scale_main(params_dict, REGION_LIST):

	params_dict = params_dict.copy()
	# Determine the type of Score file.
	if params_dict.get('genome2bit') != None:
		pass # return

	bedfile = params_dict.get('bedfile') # already converted into absolute path and could be opend by pysam.TabixFile
	tokeep = params_dict.get('keep_bases')
	fasta = params_dict.get('genomefasta')
	ref_fa = Fasta(fasta)

	#######
	# create tmp bed file to save subtracted reads
	# these tmp files would be deleted after this procedure completed.
	#######
	tmp_bedfile = str(bedfile) + '.' + get_random_string(8)
	tmpfile_list = filterBed(bedfile, ref_fa, tmp_bedfile, tokeep)
	ref_fa.close()
	#################### Done


	upstream = params_dict.get('upstream')
	downstream = params_dict.get('downstream')
	binsize = params_dict.get('binsize')
	bodylength = params_dict.get('bodylength')
	cores = params_dict.get('cores')
	totalCounts = params_dict.get('totalCounts')

	normlizeMethod = params_dict.get('norm_m')
	r_l = params_dict.get('readLength')
	g_s = params_dict.get('genomeSize')


	up_bins = upstream // binsize
	mid_bins = bodylength // binsize
	down_bins = downstream // binsize
	nbins = int(up_bins + mid_bins + down_bins)
	params_dict['nBins'] = nbins

	df_ts = pd.DataFrame(0, index=[], columns=range(nbins), dtype=float)
	df_nts = pd.DataFrame(0, index=[], columns=range(nbins), dtype=float)

	# chunk REGION_LIST into small set of REGIONS respect to given cores
	nRegion = len(REGION_LIST)
	if nRegion <= cores:
		cores = nRegion
		nCut = 1
	else:
		nCut = math.ceil(nRegion / cores)


	JOBS = []
	i = 0
	while i < nRegion:
		
		sub_REGIONS = REGION_LIST[i:i+nCut]
		sub_job = tuple([tmpfile_list[0], up_bins, mid_bins, down_bins, totalCounts, normlizeMethod, r_l, g_s, sub_REGIONS])
		JOBS.append(sub_job)
		i += nCut

	pool = multiprocessing.Pool(cores)
	res = pool.map_async(scale_run_wrapper, JOBS).get()
	pool.close()
	pool.join()

	# res : ([TS,NTS], [], [], ...)
	# TS/NTS: [[4,1,2,5,6], [5,3,2,6,7], [], ...]
	for TS, NTS in res:

		sub_TS = pd.DataFrame(np.asarray(TS, dtype=float).reshape(len(TS), nbins), columns=range(nbins))
		sub_NTS = pd.DataFrame(np.asarray(NTS, dtype=float).reshape(len(NTS), nbins), columns=range(nbins))

		# ensure the columns are same in length 
		#assert sub_TS.shape[1] == sub_NTS.shape[1] == nbins

		df_ts = pd.concat([df_ts,sub_TS], axis=0)
		df_nts = pd.concat([df_nts,sub_NTS], axis=0)

	# merge TS/NTS dataframe
	DF = pd.concat([df_ts, df_nts], axis=1)
	dim = DF.shape
	DF.index = range(dim[0])
	DF.columns = range(dim[1])

	# remove temp files
	#############################
	# Annotate for Testing
	for f in tmpfile_list:
		os.remove(f)
	#############################
	return DF, params_dict.copy()




def main(args=None):
	
	# load params and choose which mode_wrapper to run
	info('Loading data...', 'out')
	regionfile = transferPath(args.region)
	mode = args.mode
	bedfile = args.bed # mutual exclusive
	genome2bit = args.genome2bit # mutual exclusive

	genomefasta = transferPath(args.genomefasta) # genome.fasta
	keep_bases = args.keep_bases # G
	upstream = args.upstream # 2000
	downstream = args.downstream # 2000
	binsize = args.binsize # 50
	rp = str(args.reference_point).upper() # tss
	stepsize = args.stepsize # none
	bodylength = args.bodylength # none
	norm_m = args.method # RPKM
	readLength = args.readLength # none
	genomeSize = args.genomeSize # none

	# create dict to save shared params
	passing_dict = {}
	# save shared params
	passing_dict['genomefasta'] = genomefasta
	passing_dict['regionfile'] = regionfile
	passing_dict['keep_bases'] = keep_bases
	passing_dict['upstream'] = upstream
	passing_dict['downstream'] = downstream
	passing_dict['binsize'] = binsize
	passing_dict['cores'] = args.cores
	passing_dict['output'] = transferPath(args.output)
	passing_dict['norm_m'] = norm_m

	####################
	# if not output.endswith('.gz')
	# sys.exit()
	######################


	# set seperated params to None
	passing_dict['readLength'] = None
	passing_dict['genomeSize'] = None
	passing_dict['rp'] = None
	passing_dict['stepsize'] = None
	passing_dict['bodylength'] = None
	passing_dict['genome2bit'] = genome2bit
	passing_dict['bedfile'] = bedfile


	# quickly check if bedfile is valid
	if bedfile != None:
		bedfile = transferPath(bedfile)
		total_counts = 0
		try:
			tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())
			all_records = tbx.fetch()

			for _ in all_records:
				total_counts += 1

			tbx.close()
			passing_dict['totalCounts'] = total_counts
			passing_dict['bedfile'] = bedfile

		except OSError as oe:
			info('Input Bedfile must be in bgzip format and has .tbi index file', 'error')
			sys.exit()


	if norm_m == 'DEPTH':

		try:
			assert readLength != None
		except AssertionError:
			info('Normalizing Method using DEPTH, you also need to specify --readLength', 'error')
			sys.exit()

		try:
			assert int(readLength) > 0
		except:
			info('--readLength must be positive integer', 'error')
			sys.exit()

		try:
			assert genomeSize != None
		except AssertionError:
			info('Normalizing Method using DEPTH, you also need to specify --genomeSize', 'error')
			sys.exit()

		try:
			genomeSize = int(genomeSize)
		except ValueError:
			genomeSize = fetchEffectiveGenomeSize(str(genomeSize))
			if genomeSize == None:
				info('Please input a valid parameter --genomeSize', 'error')
				sys.exit()

		passing_dict['readLength'] = readLength
		passing_dict['genomeSize'] = genomeSize


	# procedure : refpoint mode
	############################
	## 2022/4/10
	## To be modified
	## binsize must be divisble by stepsize
	## and total length must be divisble by binsize 
	############################
	if mode == 'refpoint':
		info('The mode you selected is "refpoint"', 'out')
		try:
			assert stepsize <= binsize
		except AssertionError:
			info('Parameter --stepsize cannot be greater than binsize, it has been set to the same length as binsize.', 'error')
			stepsize = binsize
		except TypeError:
			stepsize = binsize

		try:
			assert binsize % stepsize == 0
		except:
			info('BinSize {:.0f} is not divisble by StepSize {:.0f}'.format(binsize, stepsize), 'error')
			sys.exit()

		passing_dict['rp'] = rp
		passing_dict['stepsize'] = stepsize

		REGION_LIST, filter_info = process_regions(passing_dict, 'refpoint')
		info(filter_info)
		# REGION_LIST = [ ['chr1',10,'+'], ...]

		# run refpoint main procedure
		info('Start run refpoint precedure ...', 'out')
		DF, res_params_Dict = refpoint_main(passing_dict, REGION_LIST)
		info('Saving results...', 'out')

		# write out
		headerLine = {
		'Mode':'refpoint',
		'ReferencePoint':res_params_Dict.get('rp'),
		'Bedfile':res_params_Dict.get('bedfile'),
		'InputTotalReads':res_params_Dict.get('totalCounts'),
		'Regionfile':res_params_Dict.get('regionfile'),
		'GenomeFasta':res_params_Dict.get('genomefasta'),
		'Genome2bit':res_params_Dict.get('genome2bit'),
		'keep_bases':res_params_Dict.get('keep_bases'),
		'Upstream':res_params_Dict.get('upstream'),
		'Downstream':res_params_Dict.get('downstream'),
		'StepSize':res_params_Dict.get('stepsize'),
		'BinSize':res_params_Dict.get('binsize'),
		'SeparatorColumn':res_params_Dict.get('nBins'),
		'NormalizingMethod':res_params_Dict.get('norm_m'),
		'ReadLength':res_params_Dict.get('readLength'),
		'GenomeSize':res_params_Dict.get('genomeSize'),
		'OutputFile':res_params_Dict.get('output'),
		}

		writeHeader(res_params_Dict.get('output'), headerLine)

		DF.to_csv(res_params_Dict.get('output'), sep='\t', header=False, index=False, mode='ab', compression='gzip', encoding='utf-8')



	# procedure : scale mode
	if mode == 'scale':
		info('The mode you selected is "scale"', 'out')

		try:
			assert bodylength != None
		except AssertionError:
			info('Computing Mode using SCALE, you also need to specify --bodylength', 'error')
			sys.exit()

		passing_dict['bodylength'] = bodylength

		REGION_LIST, filter_info = process_regions(passing_dict, 'scale')
		info(filter_info)

		# run scale main procedure
		info('Start run scale precedure ...', 'out')
		DF, res_params_Dict = scale_main(passing_dict, REGION_LIST)
		info('Saving results...', 'out')

		# write out
		headerLine = {
		'Mode':'scale',
		'Bedfile':res_params_Dict.get('bedfile'),
		'InputTotalReads':res_params_Dict.get('totalCounts'),
		'Regionfile':res_params_Dict.get('regionfile'),
		'GenomeFasta':res_params_Dict.get('genomefasta'),
		'Genome2bit':res_params_Dict.get('genome2bit'),
		'keep_bases':res_params_Dict.get('keep_bases'),
		'Upstream':res_params_Dict.get('upstream'),
		'Downstream':res_params_Dict.get('downstream'),
		'BodyLength':res_params_Dict.get('bodylength'),
		'BinSize':res_params_Dict.get('binsize'),
		'SeparatorColumn':res_params_Dict.get('nBins'),
		'NormalizingMethod':res_params_Dict.get('norm_m'),
		'ReadLength':res_params_Dict.get('readLength'),
		'GenomeSize':res_params_Dict.get('genomeSize'),
		'OutputFile':res_params_Dict.get('output'),
		}

		writeHeader(res_params_Dict.get('output'), headerLine)

		DF.to_csv(res_params_Dict.get('output'), sep='\t', header=False, index=False, mode='ab', compression='gzip', encoding='utf-8')


	info('{} Done.'.format(os.path.basename(bedfile)), 'out')

if __name__ == '__main__':
	main()
