#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import pysam, gzip
from pyfaidx import Fasta
import py2bit

from ddtools.utilities import *


def _process_region(tbx, tbit, line, colors_dict):
	
	# line = ['chr1',10,20,. ,. , +, ...]
	res = line[:]
	chrom = line[0]
	start = int(line[1])
	end = int(line[2])
	feature_type = line[3]
	#strand = line[5]
	try:
		color = line[8]
	except:
		color = None

	# fetch reads
	p, m = bedCountsinRegion(tbx, chrom, start, end)

	# fetch base context
	base_A = get_base_content(chrom, start, end, tbit, fraction=False, base='A')
	base_C = get_base_content(chrom, start, end, tbit, fraction=False, base='C')
	base_G = get_base_content(chrom, start, end, tbit, fraction=False, base='G')
	base_T = get_base_content(chrom, start, end, tbit, fraction=False, base='T')

	# add color
	col_res = colors_dict.copy()
	#if feature_type not in col_res.keys() and color != None:
	if feature_type not in col_res.keys():
		col_res[feature_type] = color

	res.extend([
		abs(start-end),
		p,
		m,
		base_A,
		base_C,
		base_G,
		base_T]
		)

	return res, col_res


def main(args=None):

	info('Loading data...', 'out')
	# check input bed
	bedfile = args.bed
	if bedfile != None:
		bedfile = transferPath(bedfile)
		total_counts = 0
		try:
			tbx = pysam.TabixFile(bedfile, mode='r', parser=pysam.asBed())
			all_records = tbx.fetch()

			for _ in all_records:
				total_counts += 1

			tbx.close()

		except OSError as oe:
			info(oe, 'error')

	# quickly check input regions
	anno_file = transferPath(args.anno_file)
	if anno_file.endswith('.gz') or anno_file.endswith('.bgz'):
		fi = gzip.open(anno_file)
		subtract_test = [fi.readline().decode().strip() for _ in range(1000)]
		fi.close()
	else:
		fi = open(anno_file)
		subtract_test = [fi.readline().strip() for _ in range(1000)]
		fi.close()

	subtract_test = [x.split('\t') for x in subtract_test if len(x) != 0 and not x.startswith('#')]
	try:
		assert len(set([len(x) for x in subtract_test])) == 1
		anno_column = len(subtract_test[0])
	except AssertionError:
		info()
		sys.exit()


	info('Generate temp file ...', 'out')
	# generate tmp bed file
	ref_fa = Fasta(transferPath(args.genomefasta))
	tbit = py2bit.open(transferPath(args.genome2bit))
	keep_bases = args.keep_bases

	tmp_bedfile = str(bedfile) + '.' + get_random_string(8)
	tmpfile_list = filterBed(bedfile, ref_fa, tmp_bedfile, keep_bases)
	ref_fa.close()

	tbx = pysam.TabixFile(tmpfile_list[0], mode='r', parser=pysam.asBed())

	#
	outfile = transferPath(args.output)
	fo = open(outfile, 'w')
	colors_dict = {}

	info('Writing ...', 'out')

	if anno_file.endswith('.gz') or anno_file.endswith('.bgz'):

		with gzip.open(anno_file) as fi:

			for line in fi:
				if line.decode().startswith('#'):
					continue

				line = line.decode().strip().split('\t')
				
				ol, colors_dict = _process_region(tbx, tbit, line, colors_dict)
				fo.write(
					'\t'.join(map(str, ol)) + '\n'
					)

	else:
		with open(anno_file) as fi:

			for line in fi:
				if line.startswith('#'):
					continue

				line = line.strip().split('\t')
				
				ol, colors_dict = _process_region(tbx, tbit, line, colors_dict)
				fo.write(
					'\t'.join(map(str, ol)) + '\n'
					)

	fo.close()
	tbx.close()
	tbit.close()


	# adding header line
	headerLine = {
		'Bedfile':str(bedfile),
		'Regionfile':str(anno_file),
		'GenomeFasta':str(transferPath(args.genomefasta)),
		'Genome2bit':str(transferPath(args.genome2bit)),
		'keep_bases':keep_bases,
		'AnnoColumn':anno_column,
		'ColorsDict':colors_dict,
		'OutputFile':str(outfile)
		}

	writeHeader(outfile, headerLine)

	# remove tmp bedfile
	for f in tmpfile_list:
		os.remove(f)

	info('{} Done.'.format(os.path.basename(bedfile)), 'out')


if __name__ == '__main__':
	main()