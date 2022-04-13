#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ddtools.utilities import *



def bin_by(x, y, nbins=10):

	bins = np.linspace(0, 1, nbins + 1)
	bins[-1] += 1
	indices = np.digitize(y, bins)
	output = []
	for i in range(1, len(bins)):
		output.append(x[indices == i])
	bins = bins[:-1]

	return output, bins


def plot_main(dataGC, dataG, smLabels, ylabel, fig_title, outFig, xtype='ALL', pt='line'):

	if xtype == 'ALL':

		fig = plt.figure(figsize=(10, 5))
		ax1 = fig.add_subplot(121)
		ax2 = fig.add_subplot(122)

		GC_y, G_y = 0, 0
		for idx, reads_per_gc in enumerate(dataGC):

			reads, GC = reads_per_gc.T
			signal_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
			to_keep = [idx for idx, x in enumerate(bin_labels) if 0.2 <= x <= 0.7]
			signal_per_gc = [signal_per_gc[x] for x in to_keep]
			bin_labels = [bin_labels[x] for x in to_keep]

			y,low,up = [],[],[]
			for i in range(len(signal_per_gc)):

				if pt == 'line':
					std_value = 0

				if pt == 'std':
					std_value = np.std(signal_per_gc[i])
					
				elif pt == 'se':
					std_value = np.std(signal_per_gc[i]) / np.sqrt(len(signal_per_gc[i]))

				meanValue = np.mean(signal_per_gc[i])
				y.append(meanValue)
				low.append(meanValue-std_value)
				up.append(meanValue+std_value)
			
			GC_y = len(y)
			x = range(1, len(y)+1)
			ax1.plot(x, y, ls='-', lw=1.5, color=DefaultColors(idx), zorder=3)
			ax1.fill_between(x, low, up, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)

		for idx, reads_per_g in enumerate(dataG):

			reads, G = reads_per_g.T
			signal_per_gc, bin_labels = bin_by(reads, G, nbins=100)
			to_keep = [idx for idx, x in enumerate(bin_labels) if 0.1 <= round(x,3) <= 0.35]
			signal_per_gc = [signal_per_gc[x] for x in to_keep]
			bin_labels = [bin_labels[x] for x in to_keep]

			y,low,up = [],[],[]
			for i in range(len(signal_per_gc)):

				if pt == 'line':
					std_value = 0

				if pt == 'std':
					std_value = np.std(signal_per_gc[i])
					
				elif pt == 'se':
					std_value = np.std(signal_per_gc[i]) / np.sqrt(len(signal_per_gc[i]))

				meanValue = np.mean(signal_per_gc[i])
				y.append(meanValue)
				low.append(meanValue-std_value)
				up.append(meanValue+std_value)

			G_y = len(y)
			x = range(1, len(y)+1)
			ax2.plot(x, y, ls='-', lw=1.5, color=DefaultColors(idx), label=smLabels[idx], zorder=3)
			ax2.fill_between(x, low, up, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)


		# adjust ax1
		ax1.set_xlabel('GC ratio',fontdict=LabelFont('label'))
		ax1.set_ylabel(ylabel, fontdict=LabelFont('label'))
		ax1.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax1_xticks = list(range(0, GC_y+1,10))
		ax1.set_xticks(ax1_xticks)
		ax1.set_xlim(0, ax1_xticks[-1]+1)
		ax1.set_xticklabels([str(x*100)+'%' for x in [0.2,0.3,0.4,0.5,0.6,0.7]], fontdict=LabelFont('tick'))

		# adjust ax2
		ax2.set_xlabel('G ratio',fontdict=LabelFont('label'))
		ax2.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax2_xticks = list(range(0, G_y+1,5))
		ax2.set_xlim(0, ax2_xticks[-1]+1)
		ax2.set_xticks(ax2_xticks)
		ax2.set_xticklabels([str(x*100)+'%' for x in [0.1,0.15,0.2,0.25,0.3,0.35]], fontdict=LabelFont('tick'))
		ax2.legend(loc='center left', bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))

		fig.suptitle(fig_title, fontdict=LabelFont('title'))
		fig.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

		return

	elif xtype == 'GC':

		fig, ax = plt.subplots(figsize=(6, 5))
		GC_y = 0

		for idx, reads_per_gc in enumerate(dataGC):

			reads, GC = reads_per_gc.T
			signal_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
			to_keep = [idx for idx, x in enumerate(bin_labels) if 0.2 <= x <= 0.7]
			signal_per_gc = [signal_per_gc[x] for x in to_keep]
			bin_labels = [bin_labels[x] for x in to_keep]

			y,low,up = [],[],[]
			for i in range(len(signal_per_gc)):

				if pt == 'line':
					std_value = 0

				if pt == 'std':
					std_value = np.std(signal_per_gc[i])
					
				elif pt == 'se':
					std_value = np.std(signal_per_gc[i]) / np.sqrt(len(signal_per_gc[i]))

				meanValue = np.mean(signal_per_gc[i])
				y.append(meanValue)
				low.append(meanValue-std_value)
				up.append(meanValue+std_value)

			GC_y = len(y)
			x = range(1, len(y)+1)
			ax.plot(x, y, ls='-', lw=1.5, color=DefaultColors(idx), label=smLabels[idx], zorder=3)
			ax.fill_between(x, low, up, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)

		ax.set_xlabel('GC ratio',fontdict=LabelFont('label'))
		ax.set_ylabel(ylabel, fontdict=LabelFont('label'))
		ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax_xticks = list(range(0,GC_y+1,10))
		ax.set_xticks(ax_xticks)
		ax.set_xlim(0, ax_xticks[-1]+1)
		ax.set_xticklabels([str(x*100)+'%' for x in [0.2,0.3,0.4,0.5,0.6,0.7]], fontdict=LabelFont('tick'))
		ax.legend(loc='center left', bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))

		fig.suptitle(fig_title, fontdict=LabelFont('title'))
		fig.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

		return


	elif xtype == 'G':

		fig, ax = plt.subplots(figsize=(6, 5))
		G_y = 0

		for idx, reads_per_g in enumerate(dataG):

			reads, G = reads_per_g.T
			signal_per_gc, bin_labels = bin_by(reads, G, nbins=100)
			to_keep = [idx for idx, x in enumerate(bin_labels) if 0.1 <= round(x,3) <= 0.35]
			signal_per_gc = [signal_per_gc[x] for x in to_keep]
			bin_labels = [bin_labels[x] for x in to_keep]

			y,low,up = [],[],[]
			for i in range(len(signal_per_gc)):

				if pt == 'line':
					std_value = 0

				if pt == 'std':
					std_value = np.std(signal_per_gc[i])
					
				elif pt == 'se':
					std_value = np.std(signal_per_gc[i]) / np.sqrt(len(signal_per_gc[i]))

				meanValue = np.mean(signal_per_gc[i])
				y.append(meanValue)
				low.append(meanValue-std_value)
				up.append(meanValue+std_value)

			G_y = len(y)
			x = range(1, len(y)+1)
			ax.plot(x, y, ls='-', lw=1.5, color=DefaultColors(idx), label=smLabels[idx], zorder=3)
			ax.fill_between(x, low, up, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)

		ax.set_xlabel('G ratio',fontdict=LabelFont('label'))
		ax.set_ylabel(ylabel, fontdict=LabelFont('label'))
		ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax_xticks = list(range(0,G_y+1,5))
		ax.set_xticks(ax_xticks)
		ax.set_xlim(0, ax_xticks[-1]+1)
		ax.set_xticklabels([str(x*100)+'%' for x in [0.1,0.15,0.2,0.25,0.3,0.35]], fontdict=LabelFont('tick'))
		ax.legend(loc='center left', bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))

		fig.suptitle(fig_title, fontdict=LabelFont('title'))
		fig.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

		return


def process_woc(df, index, regionsize):

	# input: pre-filter pd.DataFrame
	# and index to select specific sample
	# subset contains 4 columns: AP AM G C
	subset = df.iloc[:,[*index,-2,-1]]
	subset.columns = ['AP', 'AM', 'GCOUNTS', 'CCOUNTS']

	# second-filter
	# remove 
	subset = subset[subset.iloc[:,:2].sum(axis=1)>0]

	subset.GCOUNTS /= regionsize
	subset.CCOUNTS /= regionsize
	subset['GCratio'] = subset.GCOUNTS + subset.CCOUNTS
	subset['ALLCOUNTS'] = subset.iloc[:,0] + subset.iloc[:,1]

	subset_total = subset.ALLCOUNTS.sum()

	subset['CPM'] = subset.ALLCOUNTS / subset_total * 1e6

	# CPM & GCratio --> GC
	GCscores = [[x[0], x[1]] for x in zip(subset.CPM.tolist(), subset.GCratio.tolist())]

	# for G
	## AP -- GCOUNTS ; AM -- CCOUNTS
	Gscores = [[x[0]/subset_total*1e6, x[1]] for x in zip(subset.AP.tolist(), subset.GCOUNTS.tolist()) if x[0]!=0] # for plus strand
	mscores = [[x[0]/subset_total*1e6, x[1]] for x in zip(subset.AM.tolist(), subset.CCOUNTS.tolist()) if x[0]!=0] # for minus strand
	Gscores.extend(mscores)

	GCscores = np.asarray(GCscores, dtype=float)
	Gscores = np.asarray(Gscores, dtype=float)
	return GCscores, Gscores



def process_wc(df, index, regionsize):

	# input: pre-filter pd.DataFrame
	# and index to select specific sample
	# subset contains 6 columns: AP AM CP CM G C
	subset = df.iloc[:,[*index,-4,-3,-2,-1]]
	subset.columns = ['AP', 'AM', 'CP','CM', 'GCOUNTS', 'CCOUNTS']

	subset['Treatcounts'] = subset.AP + subset.AM
	subset['Controlcounts'] = subset.CP + subset.CM

	# second-filter
	# remove 
	subset = subset[(subset.Treatcounts>0) & (subset.Controlcounts>0)]
	Treatcounts_all = subset.Treatcounts.sum()
	Controlcounts_all = subset.Controlcounts.sum()
	ratio = Treatcounts_all / Controlcounts_all

	subset.GCOUNTS /= regionsize
	subset.CCOUNTS /= regionsize
	subset['GCratio'] = subset.GCOUNTS + subset.CCOUNTS

	# for GC
	subset['NormV'] = subset.Treatcounts / subset.Controlcounts / ratio
	GCscores = [[x[0], x[1]] for x in zip(subset.NormV.tolist(), subset.GCratio.tolist())]

	# for G
	Gscores = [[np.log2(a/b/ratio), c] for a,b,c in zip(subset.AP.tolist(), subset.CP.tolist(), subset.GCOUNTS.tolist()) if a!=0 and b!=0]
	mscores = [[np.log2(a/b/ratio), c] for a,b,c in zip(subset.AM.tolist(), subset.CM.tolist(), subset.CCOUNTS.tolist()) if a!=0 and b!=0]
	Gscores.extend(mscores)

	GCscores = np.asarray(GCscores, dtype=float)
	Gscores = np.asarray(Gscores, dtype=float)
	return GCscores, Gscores


def main(args=None):

	ScoreFile = transferPath(args.input)
	outFig = transferPath(args.output)
	xtype = args.xtype.upper()
	file_index = args.file_index
	sampleLabels = args.sampleLabels
	fig_title = args.title
	pt = args.plotType
	if not fig_title:
		fig_title = ''

	header = readHeader(ScoreFile)
	bedfiles = header.get('Bedfiles')
	bedlabels = [os.path.basename(x).split('.')[0] for x in bedfiles]
	control = header.get('ControlFile')
	nSamples = len(bedfiles)
	regionsize = header.get('RegionSize')

	#
	try:
		assert xtype in ['ALL','G','GC']
	except AssertionError:
		info()
		sys.exit()

	#
	fIndices = [[(i-1)*2, i*2-1] for i in range(1, nSamples+1)]
	smLabels = bedlabels[:]

	if file_index != None:
		smLabels = []

		try:
			assert max(file_index) <= nSamples
		except AssertionError:
			info()
			sys.exit()

		if len(set(file_index) - set(range(1,nSamples+1))) != 0:
			info()
			sys.exit()

		fIndices = []
		for i in file_index:
			idx_loc = [(i-1)*2, i*2-1]
			fIndices.append(idx_loc)
			smLabels.append(bedlabels[i])

	#
	if sampleLabels:

		try:
			assert len(sampleLabels) == len(smLabels)
		except AssertionError:
			info()
			sys.exit()

		smLabels = sampleLabels


	# reading data
	df = pd.read_table(ScoreFile, sep='\t', comment='#', header=None, dtype=int)

	dataGC, dataG = [], []
	if not control:
		for idx in fIndices:
			fGC, fG = process_woc(df, idx, regionsize)
			dataGC.append(fGC)
			dataG.append(fG)

		ylabel = 'CPM'
		plot_main(dataGC, dataG, smLabels, ylabel, fig_title, outFig, xtype, pt)

	else:
		for idx in fIndices:
			fGC, fG = process_wc(df, idx, regionsize)
			dataGC.append(fGC)
			dataG.append(fG)

		ylabel = 'Relative Signal'
		plot_main(dataGC, dataG, smLabels, ylabel, fig_title, outFig, xtype, pt)
		
	
	# Done.


if __name__ == '__main__':
	main()







