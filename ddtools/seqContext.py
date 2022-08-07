#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import pandas as pd
import pysam
from pyfaidx import Fasta
from functools import reduce
from multiprocessing import Pool


from ddtools.utilities import *


def list_combinations(char, keep):
	'''
	deprecated

	# this function source from itertools
	# https://docs.python.org/3/library/itertools.html#itertools.combinations_with_replacement
	# return: a list of tuples
	# [('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'T'), ('C', 'C'), ('C', 'G'), ('C', 'T'), ('G', 'G'), ('G', 'T'), ('T', 'T')]
	# or : [('A', 'A', 'A'), ('A', 'A', 'C'), ('A', 'A', 'G'), ...]
	pool = tuple(iterable)
	n = len(pool)
	if not n and r:
		return
	indices = [0] * r
	yield tuple(pool[i] for i in indices)
	while True:
		for i in reversed(range(r)):
			if indices[i] != n - 1:
				break
		else:
			return
		indices[i:] = [indices[i] + 1] * (r - i)
		yield tuple(pool[i] for i in indices)
	'''
	df = pd.DataFrame(0, index=char, columns=keep)
	for i in range(df.shape[0]):
		for j in range(df.shape[1]):
			df.iloc[i, j] = df.index[i] + df.columns[j]

	return df.values.flatten()


def cat_dict(dict_a, dict_b):

	assert sorted(dict_a.keys()) == sorted(dict_b.keys())
	res = {}
	for k in sorted(dict_a.keys()):
		merged_value = dict_a.get(k) + dict_b.get(k)
		res[k] = merged_value

	return res


def filter_seqs(seqList, region_size, left_bias, keep_bases):

	res = []
	keep_bases = keep_bases
	bases_length = len(keep_bases[0])

	for seq in seqList:

		if len(seq) != region_size:
			continue

		if bases_length == 2 and seq[left_bias-1:left_bias+1] not in keep_bases:
			continue

		if bases_length == 1 and seq[left_bias:left_bias+1] not in keep_bases:
			continue

		else:
			res.append(seq)

	return res


def getSeqs_wrapper(args):

	return getSeqs(*args)


def getSeqs(contig, region_start, region_end, left_bias, right_bias, keep_bases):

	res = []
	region_size = left_bias + right_bias + 1

	tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())
	ref_fa = Fasta(fastafile)

	records = tbx.fetch(contig, region_start, region_end)

	for record in records:

		chrom = record.contig
		start = record.start
		end = record.end
		strand = record.strand

		if strand == '+':

			re_start = start - left_bias
			re_end = end + right_bias


			sequence = fetch_sequence(ref_fa, chrom, re_start, re_end)
			if not sequence:
				continue
			try:
				sequence = sequence.seq.upper()
			except ValueError:
				print('Invalid seq: {}\t{}:{}-{}'.format(sequence,chrom,re_start,re_end))
				continue

			if sequence.count('N'):
				continue

			res.append(sequence)

		elif strand == '-':

			re_start = start - right_bias
			re_end = end + left_bias

			sequence = fetch_sequence(ref_fa, chrom, re_start, re_end)
			if not sequence:
				continue

			try:
				sequence = sequence.reverse.complement.seq.upper()
			except ValueError:
				print('Invalid seq: {}\t{}:{}-{}'.format(sequence,chrom,re_start,re_end))
				continue

			if sequence.count('N'):
				continue

			res.append(sequence)

	tbx.close()
	ref_fa.close()

	res = filter_seqs(res, region_size, left_bias, keep_bases)
	return res


def cutGenome(fa, trim_contigs):

	res = []
	trim_contigs = trim_contigs

	fa_dict = {}
	for k, v in fa.items():
		if k in trim_contigs:
			continue
		fa_dict[k] = len(v[:].seq.strip())

	#total = sum(fa_dict.values())

	chunksize = 5e6

	for chrom, max_length in fa_dict.items():

		if max_length <= chunksize:
			region_list = [chrom, 0, max_length]
			res.append(region_list)
			continue

		r = np.arange(0, max_length, chunksize)
		if r[-1] < max_length:
			r = np.append(r, max_length)

		for i in range(len(r)-1):
			region_list = [chrom, r[i], r[i+1]]
			res.append(region_list)
	
	return res


def count_multi_bases(seqList, left_bias):

	sizes = len(seqList[0])
	mono_df = pd.DataFrame(0, index=range(sizes), columns=sorted(reduce(list_combinations, [['A','C','G','T']] * 1)), dtype=int)
	di_df = pd.DataFrame(0, index=range(sizes-1), columns=sorted(reduce(list_combinations, [['A','C','G','T']] * 2)), dtype=int)
	tri_dict = {}
	for tri_type in sorted(reduce(list_combinations, [['A','C','G','T']] * 3)):
		tri_dict[tri_type] = 0

	for seq in seqList:

		for j in range(sizes):
			mono_base = seq[j]
			mono_df.loc[j, mono_base] += 1 # maybe the Limiting step ?

		for j in range(sizes-1):
			di_base = seq[j:j+2]
			di_df.loc[j, di_base] += 1

		tri_base = seq[left_bias-1:left_bias+2]
		tri_dict[tri_base] += 1

	return mono_df, di_df, tri_dict
	

###############
'''
2022.3.7 adding

The previous script was stuck at function "count_multi_bases".
So split this step into two parts   

'''

def count_by_index(seq_list, index, extend):
	base_dict = {}

	for seq in seq_list:
		mbase = seq[index:index+extend]
		k = base_dict.get(mbase, 0)
		base_dict[mbase] = k + 1

	return base_dict



def count_multiple_bases(seqList, left_bias):

	sizes = len(seqList[0])
	rb = sizes - left_bias -1

	idx_range = [(-1)*x for x in reversed(range(1,left_bias+1))]
	idx_range.extend([0])
	idx_range.extend([x for x in range(1,rb+1)])

	# mono
	mono_df = pd.DataFrame(0, index=range(sizes), columns=sorted(reduce(list_combinations, [['A','C','G','T']] * 1)), dtype=int)
	for i in range(sizes):
		bd = count_by_index(seqList, i, 1)
		for k in bd.keys():
			mono_df.loc[i, k] = bd.get(k)

	mono_df.index = idx_range

	# di
	di_df = pd.DataFrame(0, index=range(sizes-1), columns=sorted(reduce(list_combinations, [['A','C','G','T']] * 2)), dtype=int)
	for i in range(sizes-1):
		bd = count_by_index(seqList, i, 2)
		for k in bd.keys():
			di_df.loc[i, k] = bd.get(k)

	# tri
	tri_dict = {}
	for tri_type in sorted(reduce(list_combinations, [['A','C','G','T']] * 3)):
		tri_dict[tri_type] = 0

	for seq in seqList:
		tri_base = seq[left_bias-1:left_bias+2]
		tri_dict[tri_base] += 1


	return mono_df, di_df, tri_dict


###############


def main(args=None):

	# args:
	#		input(bed) prefix genome_fasta left right keep_bases threads

	info('Loading data...', 'out')
	left_bias = args.left
	right_bias = args.right
	keep_bases = args.keep_bases
	#region_size = left_bias + right_bias + 1
	#bases_length = len(keep_bases[0])
	threads = args.threads

	global bedfile, fastafile
	bedfile = transferPath(args.input)
	fastafile = transferPath(args.fasta)
	

	try:
		assert max([len(x) for x in keep_bases]) <= 2
	except AssertionError:
		info('Among bases "--keep_bases {}", the maximum length of bases cannot exceed 2'.format(' '.join(map(str,a))), 'error')
		sys.exit()

	try:
		assert len(set([len(x) for x in keep_bases])) == 1
	except AssertionError:
		info('Bases specified "--keep_bases {}" are not same in length'.format(' '.join(map(str,a))), 'error')
		sys.exit()

	try:
		tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())
		contigs = tbx.contigs
		tbx.close()
	except OSError as oe:
		info(oe, 'error')


	ref_fa = Fasta(fastafile)
	assert len(set(contigs) - set(ref_fa.keys())) == 0

	#####
	# some contigs in ref_fa were not in contigs
	#
	trim_contigs = set(ref_fa.keys()) - set(contigs)

	#####
	# cut genome into chunk
	info('Cutting genome into chunks.', 'out')
	REGIONS = cutGenome(ref_fa, trim_contigs)

	JOBS = []
	for r in REGIONS:
		r.extend([left_bias, right_bias, keep_bases])
		JOBS.append(tuple(r))

	info('Start allocate jobs and run.', 'out')
	pool = Pool(threads)
	res = pool.map_async(getSeqs_wrapper, JOBS).get()
	pool.close()
	pool.join()


	info('Collecting sequences done.', 'out')

	seqMelt = [x for sub in res for x in sub]

	info('Counting bases.', 'out')
	mono_df, di_df, tri_dict = count_multiple_bases(seqMelt, left_bias)

	# prepare three files for output
	prefix_out = args.prefix
	op_dict = {}
	for t in ['mono', 'di', 'tri']:
		filename = prefix_out + '.' + t + '.txt'
		op_dict[t] = FileIO(filename, None, 'file')

	info('Saving results...', 'out')
	# save results
	## mono
	mono_header_line = {
	'type':'mono',
	'keep_bases':keep_bases,
	'left_bias':left_bias,
	'right_bias':right_bias
	}
	
	writeHeader(op_dict.get('mono'), mono_header_line)
	mono_df.to_csv(op_dict.get('mono'), sep='\t', header=True, index=True, quoting=None, mode='a')

	## di
	di_header_line = {
	'type':'di',
	'keep_bases':keep_bases
	}

	writeHeader(op_dict.get('di'), di_header_line)
	di_df.to_csv(op_dict.get('di'), sep='\t', header=True, index=True, quoting=None, mode='a')

	## tri
	tri_header_line = {
	'type':'tri',
	'keep_bases':keep_bases
	}

	writeHeader(op_dict.get('tri'), tri_header_line)
	with open(op_dict.get('tri'), 'a') as fo:
		for k in sorted(tri_dict.keys()):
			v = tri_dict.get(k)
			fo.write('{}\t{:.0f}\n'.format(k, v))


	info('{} Done.'.format(os.path.basename(bedfile)), 'out')

	ref_fa.close()


if __name__ == '__main__':
	main()



