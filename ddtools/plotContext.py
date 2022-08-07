#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ddtools.seqContext import list_combinations
from ddtools.utilities import *


def _groupByLastBase(indict):

	#######
	# params:tokeep defined by user
	# length of indict.keys() must be multiples of 4
	# for each spot, the posibilities are 4, 8, 12 or 16
	#######

	def colorbybase(i):

		BaseColors = {
			'A':'#43CD80',
			'C':'#008B00',
			'G':'#1E90FF',
			'T':'#4682B4',
		}

		return BaseColors.get(i)

	features = list(indict.keys())
	n = len(features)

	res_ticks = []
	res_ticklabels = []
	res_heights = []
	res_colors = []

	# sort base by its potential energy
	sort_rules = {'T':0,'C':1,'A':2,'G':3}
	sort_features = sorted(features, key=lambda x:(sort_rules[x[-1]], sort_rules[x[0]], sort_rules[x[1]]))

	if n == 16:
		# each group has 4 elements
		# around x-tick == 1
		# range: 0.5 - 1.5
		# step of groups: 1.5
		# width: 0.3
		width = 0.3
		base_ticks = np.linspace(0.5, 1.5, 4)
		for i in range(4):
			res_ticks.extend(list(base_ticks + i*1.5))

		for i, t in enumerate(sort_features):
			res_ticklabels.append(t)
			v = indict.get(t)
			res_heights.append(v)
			res_colors.append(colorbybase(t[-1]))

		return res_ticks, res_ticklabels, res_heights, res_colors, width


	elif n == 32:
		# each group has 8 elements
		# around x-tick == 1
		# range: 0.5 - 2.5
		# step of groups: 3
		# width: 0.24
		width = 0.24
		base_ticks = np.linspace(0.5, 2.5, 8)
		for i in range(4):
			res_ticks.extend(list(base_ticks + i*2.5))

		for i, t in enumerate(sort_features):
			res_ticklabels.append(t)
			v = indict.get(t)
			res_heights.append(v)
			res_colors.append(colorbybase(t[-1]))

		return res_ticks, res_ticklabels, res_heights, res_colors, width


	elif n == 48:
		# each group has 12 elements
		# around x-tick == 1
		# range: 0.5 - 2.5
		# step of groups: 3
		# width: 0.15
		width = 0.15
		base_ticks = np.linspace(0.5, 2.5, 12)
		for i in range(4):
			res_ticks.extend(list(base_ticks + i*2.5))

		for i, t in enumerate(sort_features):
			res_ticklabels.append(t)
			v = indict.get(t)
			res_heights.append(v)
			res_colors.append(colorbybase(t[-1]))

		return res_ticks, res_ticklabels, res_heights, res_colors, width


	elif n == 64:
		# each group has 16 elements
		# around x-tick == 1
		# range: 0.5 - 2.5
		# step of groups: 3
		# width: 0.15
		width = 0.15
		base_ticks = np.linspace(0.5, 3.5, 16)
		for i in range(4):
			res_ticks.extend(list(base_ticks + i*3.5))

		for i, t in enumerate(sort_features):
			res_ticklabels.append(t)
			v = indict.get(t)
			res_heights.append(v)
			res_colors.append(colorbybase(t[-1]))

		return res_ticks, res_ticklabels, res_heights, res_colors, width

	else:
		info()
		sys.exit()


def plot_mono(infile, outfile, figType, lb, rb, title):
	
	# loading data
	df = pd.read_table(infile, comment='#', header=0, index_col=0, dtype=int)
	x = range(lb+rb+1)
	xticklabels = [-1*i for i in range(lb, 0, -1)] + [0] + [i for i in range(1, rb+1, 1)]

	total_count = int(df.sum(axis=1)[0])
	a = np.array([x / total_count for x in df['A'].tolist()])
	c = np.array([x / total_count for x in df['C'].tolist()])
	g = np.array([x / total_count for x in df['G'].tolist()])
	t = np.array([x / total_count for x in df['T'].tolist()])

	fig, ax = plt.subplots(figsize=(6,5))

	if figType == 'bar':

		bar1 = ax.bar(x, a, color='#43CD80', label='A', alpha=1, zorder=2)
		bar2 = ax.bar(x, c, bottom=a, color='#008B00', label='C', alpha=1, zorder=2)
		bar3 = ax.bar(x, g, bottom=a + c, color='#1E90FF', label='G', alpha=1, zorder=2)
		bar4 = ax.bar(x, t, bottom=a + c + g, color='#4682B4', label='T', alpha=1, zorder=2)

		#ax.legend(loc='lower left', ncol=4, bbox_to_anchor=(0, 1.02, 1, 0.2), mode='expand')
		ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))
		ax.tick_params(bottom=False)
		ax.tick_params(left=False)
		ax.set_xticks(x)
		ax.set_xticklabels(xticklabels, fontdict=LabelFont('tick'))
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.set_ylim(0, 1)
		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		ax.set_yticks([0, .2, .4, .6, .8, 1])
		ax.set_yticklabels([str(j) + '%' for j in [int(i*100) for i in [0, .2, .4, .6, .8, 1]]], fontdict=LabelFont('tick'))
		ax.grid(axis='y', ls='-', lw=1.5, color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return

	elif figType == 'line':

		ax.plot(x, a, color='#43CD80', ls='--', lw=1.5, label='A', marker='o', ms=5, zorder=2)
		ax.plot(x, c, color='#008B00', ls='--', lw=1.5, label='C', marker='o', ms=5, zorder=2)
		ax.plot(x, g, color='#1E90FF', ls='--', lw=1.5, label='G', marker='o', ms=5, zorder=2)
		ax.plot(x, t, color='#4682B4', ls='--', lw=1.5, label='T', marker='o', ms=5, zorder=2)
		ax.legend(loc='upper right', frameon=False, prop=LabelFont('legend'))
		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		ax.set_ylim(0,1)
		ax.set_yticks([0, .2, .4, .6, .8, 1])
		ax.set_yticklabels([str(j) + '%' for j in [int(i*100) for i in [0, .2, .4, .6, .8, 1]]], fontdict=LabelFont('tick'))
		ax.set_xticks(x)
		ax.set_xticklabels(xticklabels, fontdict=LabelFont('tick'))
		ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return


def plot_di(infile, outfile, tokeep, figType, title):
	
	# loading data
	df = pd.read_table(infile, comment='#', header=0, index_col=0, dtype=int)
	total_count = int(df.sum(axis=1)[0])
	df /= total_count

	if not tokeep:
		tokeep = df.columns.tolist()
		tokeep.sort()

	try:
		assert len(set(tokeep) - set(df.columns)) == 0
	except AssertionError:
		info()
		sys.exit()

	x = range(df.shape[0])
	xticklabels = []
	for i in range(df.shape[0]):
		xticklabels.append('%s-%s' % (i+1, i+2))


	fig, ax = plt.subplots(figsize=(6,5))

	bar_base = np.array([0 for _ in x])

	if figType == 'bar':

		for i in range(len(tokeep)):
			k = tokeep[i]
			pct = np.asarray(df[k])
			ax.bar(x, pct, bottom=bar_base, label=k, color=DefaultColors(i), zorder=2)
			bar_base = bar_base + pct

		ax.tick_params(bottom=False)
		ax.tick_params(left=False)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.set_ylim(0,1)
		ax.set_yticks([0, .2, .4, .6, .8, 1])
		ax.set_yticklabels([str(j) + '%' for j in [int(i*100) for i in [0, .2, .4, .6, .8, 1]]], fontdict=LabelFont('tick'))
		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))
		ax.set_xticks(x)
		ax.set_xticklabels(xticklabels, fontdict=LabelFont('tick'))
		ax.grid(axis='y', lw=1.5, color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return

	elif figType == 'line':

		for i in range(len(tokeep)):
			k = tokeep[i]
			pct = np.asarray(df[k])
			ax.plot(x, pct, label=k, ls='--', lw=1, marker='o', color=DefaultColors(i), ms=5, zorder=2)

		ax.legend(loc='upper right', frameon=False, prop=LabelFont('legend'))
		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		ax.set_ylim(0,1)
		ax.set_yticks([0, .2, .4, .6, .8, 1])
		ax.set_yticklabels([str(j) + '%' for j in [int(i*100) for i in [0, .2, .4, .6, .8, 1]]], fontdict=LabelFont('tick'))
		ax.set_xticks(x)
		ax.set_xticklabels(xticklabels, fontdict=LabelFont('tick'))
		ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return


def plot_tri(infile, outfile, tokeep, title, cfile=None):

	triDict = _filereading(infile)
	keys = list(triDict.keys())
	keys.sort()
	fkeys = [x for x in keys if x[1] in tokeep]
	t_total = sum([triDict.get(i) for i in fkeys])

	if cfile:
		cDict = _filereading(cfile)
		cfkeys = [x for x in cDict.keys() if x[1] in tokeep]
		c_total = sum([cDict.get(i) for i in cfkeys])

		try:
			assert len(set(fkeys) - set(cDict.keys())) == 0
		except AssertionError:
			info()
			sys.exit()

		# preprocessing
		t2c = {}
		for key in fkeys:
			try:
				_ratio = (triDict.get(key) * c_total) / (cDict.get(key) * t_total)
			except ZeroDivisionError:
				info()
				sys.exit()

			t2c[key] = np.log2(_ratio)

		#
		res_ticks, res_ticklabels, res_heights, res_colors, res_width = _groupByLastBase(t2c)

		fig, ax = plt.subplots(figsize=(6,5))
		ax.bar(res_ticks, res_heights, color=res_colors,width=res_width, zorder=2)

		ax.tick_params(bottom=False)
		ax.tick_params(left=False)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.spines['left'].set_color('none')

		ax.set_xticks(res_ticks)
		ax.set_xticklabels(res_ticklabels, fontdict=LabelFont('tick'), rotation=90)

		fig_yticks = list(ax.get_yticks())
		ax.set_ylim(fig_yticks[0], fig_yticks[-1])
		ax.set_yticks(fig_yticks)
		ax.set_yticklabels(fig_yticks, fontdict=LabelFont('tick'))

		ax.grid(axis='y', lw=1, color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)
		ax.set_ylabel('Relative Signal', fontdict=LabelFont('label'))

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return

	else:

		res_heights = []
		x = range(len(fkeys))
		for k in fkeys:
			res_heights.append(triDict.get(k))


		res_heights = np.asarray(res_heights, dtype=float) / t_total

		fig, ax = plt.subplots(figsize=(6,5))
		ax.bar(x, res_heights, zorder=2)

		ax.tick_params(bottom=False)
		ax.tick_params(left=False)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.spines['left'].set_color('none')

		ax.set_xticks(x)
		ax.set_xticklabels(fkeys, fontdict=LabelFont('tick'), rotation=90)

		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		fig_yticks = list(ax.get_yticks())
		fig_yticks_limits = [fig_yticks[0], fig_yticks[-1]*1.2]
		re_yticks = np.linspace(fig_yticks_limits[0], fig_yticks_limits[1], 5)
		re_yticklabels = [str(j) + '%' for j in [int(i*100) for i in re_yticks]]
		
		ax.set_ylim(fig_yticks_limits[0], fig_yticks_limits[-1])
		ax.set_yticks(re_yticks)
		ax.set_yticklabels(re_yticklabels, fontdict=LabelFont('tick'))
		ax.grid(axis='y', lw=1, color='grey', zorder=1)
		ax.set_title(title, fontdict=LabelFont('title'), pad=10)

		plt.savefig(outfile, dpi=300, bbox_inches='tight')
		plt.close()

		return


def _filereading(f):

	res = {}

	with open(f) as fi:

		for line in fi:
			if line.startswith('#'):
				continue

			if len(line.strip()) == 0:
				continue

			else:
				line = line.strip().split('\t')
				res[line[0]] = int(line[1])

	return res



def main(args=None):

	# load params
	seqFile = transferPath(args.input)
	outFig = transferPath(args.output)
	figTitle = args.title
	if not figTitle:
		figTitle = os.path.basename(seqFile).split('.')[0]
	# parse args

	# parse header line

	'''
	mono_header_line = {
	'type':'mono',
	'keep_bases':keep_bases,
	'left_bias':left_bias,
	'right_bias':right_bias
	}

	## di
	di_header_line = {
	'type':'di',
	'keep_bases':keep_bases
	}

	## tri
	tri_header_line = {
	'type':'tri',
	'keep_bases':keep_bases
	}
	'''



	##############
	# 2022/5/3
	# to be modified
	# 检查此处读入文件的keepbases和参数keepbase是否包含关系
	###############

	params_dict = readHeader(seqFile)
	seqType = params_dict.get('type')
	#

	if seqType == 'mono':
		left_bias = params_dict.get('left_bias')
		right_bias = params_dict.get('right_bias')
		figType = args.mono_type

		plot_mono(seqFile, outFig, figType, left_bias, right_bias, figTitle)

	if seqType == 'di':
		figType = args.di_type
		di_keep = args.di_keep

		plot_di(seqFile, outFig, di_keep, figType, figTitle)

	if seqType == 'tri':
		tri_keep = args.tri_keep
		controlFile = None
		if args.control:
			controlFile = transferPath(args.control)

		plot_tri(seqFile, outFig, tri_keep, figTitle, controlFile)


if __name__ == '__main__':
	main()









		






				






	
