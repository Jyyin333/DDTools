#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ddtools.utilities import *
from ddtools.seqContext import list_combinations


def plot_main(data, colors, backbase, outFig, figTitle, sort_keys):

	# treat only
	# barplot
	if data.shape[1] == 7:

		x, y = [], []
		dataColors = []
		bgb_mean, bgb_yerr = [], []

		for idx, ft in enumerate(sort_keys):
			x.append(idx+1)
			subset = data[data.Feature_type == ft]
			ft_counts = subset.gTreat.sum()
			y.append(ft_counts)
			dataColors.append(colors.get(ft))

			if backbase and len(backbase) == 2:
				ft_backbase_ratio = subset.loc[:][backbase].sum(axis=1) / subset.Length
				ft_backbase_mean = ft_backbase_ratio.mean()
				ft_backbase_se = np.std(ft_backbase_ratio) / np.sqrt(len(ft_backbase_ratio))
				bgb_mean.append(ft_backbase_mean)
				bgb_yerr.append(ft_backbase_se)

			if backbase and len(backbase) == 1:
				ft_backbase_ratio = subset.loc[:][backbase[0]] / subset.Length
				ft_backbase_mean = ft_backbase_ratio.mean()
				ft_backbase_se = np.std(ft_backbase_ratio) / np.sqrt(len(ft_backbase_ratio))
				bgb_mean.append(ft_backbase_mean)
				bgb_yerr.append(ft_backbase_se)

		######################
		# 2022/4/11
		# To be modified
		# Divide the values in Y by the proportions of each region in case of one type of region-proportion is too large
		######################


		allftcounts = sum(y)
		y = [x / allftcounts for x in y]

		fig, ax = plt.subplots(figsize=(6,5))
		ax.bar(x, y, color=dataColors, width=0.8, zorder=2)

		ax.tick_params(bottom=False)
		ax.tick_params(left=False)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.spines['left'].set_color('none')

		ax.set_xticks(x)
		ax.set_xticklabels(sort_keys, fontdict=LabelFont('tick'))

		ax.set_ylabel('% Percentage', fontdict=LabelFont('label'))
		fig_yticks = list(ax.get_yticks())
		fig_yticks_limits = [fig_yticks[0], fig_yticks[-1]*1.2]
		re_yticks = np.linspace(fig_yticks_limits[0], fig_yticks_limits[1], 5)
		re_yticklabels = [str(j) + '%' for j in [int(i*100) for i in re_yticks]]

		ax.set_ylim(fig_yticks_limits[0], fig_yticks_limits[-1])
		ax.set_yticks(re_yticks)
		ax.set_yticklabels(re_yticklabels, fontdict=LabelFont('tick'))
		ax.grid(axis='y', lw=1, color='grey', zorder=1)
		ax.set_title(figTitle, fontdict=LabelFont('title'), pad=10)

		if backbase:
			ax2 = ax.twinx()
			ax2.tick_params(bottom=False)
			ax2.errorbar(x, bgb_mean, xerr=0.25, yerr=bgb_yerr, fmt='o', ecolor='#778899', elinewidth=1.5, zorder=4, mfc='#B22222', mec='#B22222')
			ax2_fig_yticks = list(ax2.get_yticks())
			ax2_re_yticklabels = [str(j) + '%' for j in [int(round(i * 100, 1)) for i in ax2_fig_yticks]]
			ax2.set_yticks(ax2_fig_yticks)
			ax2.set_yticklabels(ax2_re_yticklabels, fontdict=LabelFont('tick'))
			ax2.set_ylabel('% {} PCT'.format(''.join(map(str, backbase))), fontdict=LabelFont('label'))


		plt.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

		return


	# contains control
	# boxplot
	else:
		x, y = [], []
		dataColors = []
		bgb_mean, bgb_yerr = [], []

		for idx, ft in enumerate(sort_keys):
			x.append(idx+1)
			subset = data[data.Feature_type == ft]
			y.append(subset.Log2Norm.tolist())
			dataColors.append(colors.get(ft))

			if backbase and len(backbase) == 2:
				ft_backbase_ratio = subset.loc[:][backbase].sum(axis=1) / subset.Length
				ft_backbase_mean = ft_backbase_ratio.mean()
				ft_backbase_se = np.std(ft_backbase_ratio) / np.sqrt(len(ft_backbase_ratio))
				bgb_mean.append(ft_backbase_mean)
				bgb_yerr.append(ft_backbase_se)

			if backbase and len(backbase) == 1:
				ft_backbase_ratio = subset.loc[:][backbase[0]] / subset.Length
				ft_backbase_mean = ft_backbase_ratio.mean()
				ft_backbase_se = np.std(ft_backbase_ratio) / np.sqrt(len(ft_backbase_ratio))
				bgb_mean.append(ft_backbase_mean)
				bgb_yerr.append(ft_backbase_se)


		fig, ax = plt.subplots(figsize=(6,5))
		ax.tick_params(bottom=False)
		bp = ax.boxplot(y, patch_artist=True, labels=None, notch=0, whis=0.45, showcaps=True, showfliers=False)
		plt.setp(bp['whiskers'], color='black', linestyle='dashed',linewidth=0.8)
		plt.setp(bp['medians'], color='black',linewidth=0.6)
		[bp['boxes'][i].set(facecolor=dataColors[i], alpha=0.8) for i in range(len(y))]

		ax.set_xticks(x)
		ax.set_xticklabels(sort_keys, fontdict=LabelFont('tick'))
		ax.set_title(figTitle, fontdict=LabelFont('title'), pad=10)
		ax.set_ylabel('Relative Signal', fontdict=LabelFont('label'))

		if backbase:
			ax2 = ax.twinx()
			ax2.tick_params(bottom=False)
			ax2.errorbar(x, bgb_mean, xerr=0.25, yerr=bgb_yerr, fmt='o', ecolor='#778899', elinewidth=1.5, zorder=4, mfc='#B22222', mec='#B22222')
			ax2_fig_yticks = list(ax2.get_yticks())
			ax2_re_yticklabels = [str(j) + '%' for j in [int(round(i * 100, 1)) for i in ax2_fig_yticks]]
			ax2.set_yticks(ax2_fig_yticks)
			ax2.set_yticklabels(ax2_re_yticklabels, fontdict=LabelFont('tick'))
			ax2.set_ylabel('% {} PCT'.format(''.join(map(str, backbase))), fontdict=LabelFont('label'))


		plt.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

		return


def prepareDatawoc(treatGroup):

	res = {}

	treatCounts = treatGroup.Plus_count + treatGroup.Minus_counts

	mdf = pd.DataFrame({
			'Feature_type':treatGroup.Feature_type,
			'gTreat':treatCounts,
			'Length':treatGroup.Length,
			'A':treatGroup.A,
			'C':treatGroup.C,
			'G':treatGroup.G,
			'T':treatGroup['T']})

	return mdf


def prepareDatawc(treatGroup, controlGroup):

	res = {}

	treatCounts = treatGroup.Plus_count + treatGroup.Minus_counts


	controlCounts = controlGroup.Plus_count + controlGroup.Minus_counts

	mdf = pd.DataFrame({
			'Feature_type':treatGroup.Feature_type,
			'gTreat':treatCounts,
			'gControl':controlCounts,
			'Length':treatGroup.Length,
			'A':treatGroup.A,
			'C':treatGroup.C,
			'G':treatGroup.G,
			'T':treatGroup['T']})

	# filtering
	mdf = mdf[(mdf.gTreat > 0) & (mdf.gControl > 0)]
	mdf['Norm'] = (mdf.gTreat / mdf.gControl) * (mdf.gControl.sum() / mdf.gTreat.sum())
	mdf['Log2Norm'] = np.log2(mdf.Norm)

	return mdf



def main(args=None):

	annoFile = transferPath(args.input)
	outFig = transferPath(args.output)
	userColors = args.colors
	bgb = args.background_base
	resort = args.resort
	figTitle = args.title
	if not figTitle:
		figTitle = os.path.basename(annoFile).split('.')[0]

	# check bgb validity
	if bgb:
		if len(bgb) > 2:
			info()
			sys.exit()

		if bgb not in ['A','C','G','T'] and bgb not in list_combinations(['A','C','G','T'], ['A','C','G','T']):
			info()
			sys.exit()

		bgb = [x for x in bgb]

	header = readHeader(annoFile)
	ColorsDict = header.get('ColorsDict')
	sort_keys = sorted(list(ColorsDict.keys()))
	# check resort validity
	if resort:
		try:
			assert sorted(resort) == sort_keys
			sort_keys = resort
		except AssertionError:
			info()
			sys.exit()

	if userColors:
		ckeys = list(ColorsDict.keys())
		try:
			assert len(userColors) == len(ckeys)
		except:
			info()
			sys.exit()

		ColorsDict = {}
		for i, k in enumerate(ckeys):
			ColorsDict[k] = userColors[i]
		
	adjColors = {}
	color_idx = 0
	for k, v in ColorsDict.items():

		if v == None:
			adjColors[k] = DefaultColors(color_idx)
			color_idx += 1
			continue

		first_code = v[0]

		if first_code == '#' and len(v) == 7 and v[1:].isalnum():
			# hexadecimal
			adjColors[k] = v.upper()
			continue

		if first_code.isalpha() and v.isalpha():
			# english name
			adjColors[k] = v
			continue

		if first_code.isnumeric() and ''.join(map(str,v.split(','))).isnumeric():
			# RGB value
			adjColors[k] = RGB2hexadecimal(v)
			continue


	# reading data
	endColumn = header.get('AnnoColumn')
	usecols = [3] + list(range(endColumn, endColumn+7))
	colNames = ['Feature_type', 'Length', 'Plus_count', 'Minus_counts', 'A', 'C', 'G', 'T']
	df = pd.read_table(annoFile, sep='\t', comment='#', header=None, usecols=usecols, names=colNames)

	if args.control:

		controlFile = transferPath(args.control)
		cheader = readHeader(controlFile)
		try:
			assert cheader.get('Regionfile') == header.get('Regionfile') and cheader.get('AnnoColumn') == endColumn
		except AssertionError:
			info()
			sys.exit()

		cdf = pd.read_table(controlFile, sep='\t', comment='#', header=None, usecols=usecols, names=colNames)
		try:
			assert cdf.shape[0] == df.shape[0]
		except AssertionError:
			info()
			sys.exit()

		mdf = prepareDatawc(df, cdf)

	else:
		mdf = prepareDatawoc(df)


	# plot
	plot_main(mdf, adjColors, bgb, outFig, figTitle, sort_keys)


if __name__ == '__main__':
	main()



