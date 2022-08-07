#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import copy

from ddtools.utilities import *


def _preprocessing_singlesample(mtx, ref_label, scale_label, pt='line'):

	res = {}
	ref_label = ref_label
	scale_label = scale_label

	header = readHeader(mtx)
	mtx_method = header.get('NormalizingMethod')
	mtx_filename = header.get('OutputFile')
	mtx_mode = header.get('Mode')
	if not ref_label:
		ref_label = header.get('ReferencePoint')

	res['fname'] = os.path.basename(mtx_filename)
	res['mode'] = mtx_mode
	res['y_label'] = mtx_method
	res['headerinfo'] = copy.deepcopy(header)

	####################
	# 2022/4/10
	# To be modified
	# se = np.std(ma, axis=0)/ np.sqrt(ma.shape[0])
	# std = np.std(ma, axis=0)
	# fill = summary(mean) +- std
	####################

	if mtx_mode == 'refpoint':

		# create xticks and its labels
		upstream = header.get('Upstream')
		downstream = header.get('Downstream')
		stepsize = header.get('StepSize')
		ref_loc = upstream // stepsize -1
		separator = header.get('SeparatorColumn')
		mtx_xticks = [0, ref_loc, separator-1]

		mtx_xticklabels = []
		start_label = '-' + str(upstream / 1000) + 'kb'
		end_label = '+' + str(downstream / 1000) + 'kb'
		mtx_xticklabels.extend([start_label, end_label])
		mtx_xticklabels.insert(1, ref_label)

		# reading data and calculate mean values
		mtx_df = pd.read_table(mtx, sep='\t', header=None, index_col=None, comment='#', dtype=float)
		mean_values = np.asarray(mtx_df.mean(axis=0).tolist(), dtype=float)

		# split dataframe into two parts
		# in computeMTX outputs, the first part of columns is TS by default
		#mtx_TS_df = mtx_df.iloc[:,:separator]
		#mtx_NTS_df = mtx_df.iloc[:,separator:mtx_df.shape[1]]
		TS_mean = mean_values[:separator]
		NTS_mean = mean_values[separator:]

		# 2022.4.10 adding
		std_value = None
		if pt == 'std':
			std_value = np.std(mtx_df, axis=0)
			
		elif pt == 'se':
			std_value = np.std(mtx_df, axis=0) / np.sqrt(mtx_df.shape[0])

		res['std'] = std_value
		res['xticks'] = mtx_xticks
		res['xticklabels'] = mtx_xticklabels
		res['ts_score'] = TS_mean
		res['nts_score'] = NTS_mean

		return res

	elif mtx_mode == 'scale':

		# create xticks and its labels
		upstream = header.get('Upstream')
		downstream = header.get('Downstream')
		binsize = header.get('BinSize')
		bodylength = header.get('BodyLength')
		separator = header.get('SeparatorColumn')
		mtx_xticks = [0, upstream//binsize - 1, (upstream+bodylength)//binsize - 1, separator-1]

		mtx_xticklabels = []
		start_label = '-' + str(upstream / 1000) + 'kb'
		end_label = '+' + str(downstream / 1000) + 'kb'
		mtx_xticklabels.extend([start_label, end_label])
		mtx_xticklabels.insert(1, scale_label[0])
		mtx_xticklabels.insert(2, scale_label[1])

		# reading data and calculate mean values
		mtx_df = pd.read_table(mtx, sep='\t', header=None, index_col=None, comment='#', dtype=float)
		mean_values = np.asarray(mtx_df.mean(axis=0).tolist(), dtype=float)
		# split dataframe into two parts
		# in computeMTX outputs, the first part of columns is TS by default
		#mtx_TS_df = mtx_df.iloc[:,:separator]
		#mtx_NTS_df = mtx_df.iloc[:,separator:mtx_df.shape[1]]
		TS_mean = mean_values[:separator]
		NTS_mean = mean_values[separator:]

		# 2022.4.10 adding
		std_value = None
		if pt == 'std':
			std_value = np.std(mtx_df, axis=0)
			
		elif pt == 'se':
			std_value = np.std(mtx_df, axis=0) / np.sqrt(mtx_df.shape[0])

		res['std'] = std_value
		res['xticks'] = mtx_xticks
		res['xticklabels'] = mtx_xticklabels
		res['ts_score'] = TS_mean
		res['nts_score'] = NTS_mean

		return res


def split_plot(nrow, ncol, treatGroup, controlGroup=None):

	figsize_x = ncol * 2.2
	figsize_y = nrow * 1.5
	fig = plt.figure(figsize=(2.5+figsize_x, 3+figsize_y))
	gs = gridspec.GridSpec(nrow, ncol, figure=fig)

	# set layout
	for i in range(nrow):

		for j in range(ncol):

			fig.add_subplot(gs[i, j])

	axs = fig.axes

	SampleNames = list(treatGroup.keys())
	# add data
	for idx, sm in enumerate(SampleNames):
		sm_dict = treatGroup.get(sm)
		sm_ts = sm_dict.get('ts_score')
		sm_nts = sm_dict.get('nts_score')
		sm_std = sm_dict.get('std')
		sm_label = sm_dict.get('sample_label')
		sm_xticks = sm_dict.get('xticks')
		sm_xticklabels = sm_dict.get('xticklabels')
		sm_ylabel = sm_dict.get('y_label')

		# if control exist, normalize to control
		if controlGroup:
			ctl_dict = controlGroup.get(sm)
			c_ts = ctl_dict.get('ts_score')
			c_nts = ctl_dict.get('nts_score')
			sm_ylabel = 'Relative Signal'

			sm_ts /= c_ts
			sm_nts /= c_nts

		x = range(len(sm_ts))
		ax = axs[idx]
		ax.plot(x, sm_ts, label='TS', ls='-', lw=1.5, color='#0000FF', zorder=3)
		ax.plot(x, sm_nts, label='NTS', ls='-', lw=1.5, color='#FF0000', zorder=3)

		if sm_std:
			sm_ts_lower = sm_ts-sm_std
			sm_ts_upper = sm_ts+sm_std
			sm_nts_lower = sm_nts-sm_std
			sm_nts_upper = sm_nts+sm_std

			if controlGroup:
				ctl_dict = controlGroup.get(sm)
				c_ts = ctl_dict.get('ts_score')
				c_nts = ctl_dict.get('nts_score')
				c_std = ctl_dict.get('std')
				c_ts_lower = c_ts-c_std
				c_ts_upper = c_ts+c_std
				c_nts_lower = c_nts-c_std
				c_nts_upper = c_nts+c_std

				sm_ts_lower /= c_ts_lower
				sm_ts_upper /= c_ts_upper
				sm_nts_lower /= c_nts_lower
				sm_nts_upper /= c_nts_upper

			ax.fill_between(x, sm_ts_lower, sm_ts_upper, facecolor='#0000FF', alpha=0.5, zorder=2)
			ax.fill_between(x, sm_nts_lower, sm_nts_upper, facecolor='#FF0000', alpha=0.5, zorder=2)

		ax.set_xticks(sm_xticks)
		ax.set_xticklabels(sm_xticklabels, fontdict=LabelFont('tick'))
		ax.set_xlabel(sm_label, fontdict=LabelFont('label'))
		ax.set_ylabel(sm_ylabel, fontdict=LabelFont('label'))
		ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
		ax.legend(loc='center left', bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))

	return fig




def merge_plot(treatGroup, controlGroup=None):

	fig = plt.figure(figsize=(6, 5))
	gs = gridspec.GridSpec(1, 2, figure=fig)
	
	# set layout
	for i in range(1):

		for j in range(2):

			fig.add_subplot(gs[i, j])

	axs = fig.axes
	# in this context
	# layout of figure must contains one row and two columns
	ts_ax = axs[0]
	nts_ax = axs[1]

	ts_merge, nts_merge = [], []
	SampleNames = list(treatGroup.keys())

	# add data
	for idx, sm in enumerate(SampleNames):
		sm_dict = treatGroup.get(sm)
		sm_ts = sm_dict.get('ts_score')
		sm_nts = sm_dict.get('nts_score')
		sm_std = sm_dict.get('std')
		sm_label = sm_dict.get('sample_label')
		sm_xticks = sm_dict.get('xticks')
		sm_xticklabels = sm_dict.get('xticklabels')
		sm_ylabel = sm_dict.get('y_label')

		# if control exist, normalize to control
		if controlGroup:
			ctl_dict = controlGroup.get(sm)
			c_ts = ctl_dict.get('ts_score')
			c_nts = ctl_dict.get('nts_score')
			sm_ylabel = 'Relative Signal'

			sm_ts /= c_ts
			sm_nts /= c_nts


		x = range(len(sm_ts))
		ts_ax.plot(x, sm_ts, ls='-', lw=1.2, color=DefaultColors(idx), zorder=3)
		nts_ax.plot(x, sm_nts, ls='-', lw=1.2, color=DefaultColors(idx), zorder=3, label=sm_label)

		if sm_std:
			sm_ts_lower = sm_ts-sm_std
			sm_ts_upper = sm_ts+sm_std
			sm_nts_lower = sm_nts-sm_std
			sm_nts_upper = sm_nts+sm_std

			if controlGroup:
				ctl_dict = controlGroup.get(sm)
				c_ts = ctl_dict.get('ts_score')
				c_nts = ctl_dict.get('nts_score')
				c_std = ctl_dict.get('std')
				c_ts_lower = c_ts-c_std
				c_ts_upper = c_ts+c_std
				c_nts_lower = c_nts-c_std
				c_nts_upper = c_nts+c_std

				sm_ts_lower /= c_ts_lower
				sm_ts_upper /= c_ts_upper
				sm_nts_lower /= c_nts_lower
				sm_nts_upper /= c_nts_upper

			ts_ax.fill_between(x, sm_ts_lower, sm_ts_upper, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)
			nts_ax.fill_between(x, sm_nts_lower, sm_nts_upper, facecolor=DefaultColors(idx), alpha=0.5, zorder=2)

		if idx == len(SampleNames) - 1:
			ts_ax.set_ylabel(sm_ylabel, fontdict=LabelFont('label'))
			ts_ax.set_xticks(sm_xticks)
			ts_ax.set_xticklabels(sm_xticklabels, fontdict=LabelFont('tick'))
			nts_ax.set_xticks(sm_xticks)
			nts_ax.set_xticklabels(sm_xticklabels, fontdict=LabelFont('tick'))

	ts_ax.set_xlabel('TS', fontdict=LabelFont('label'))
	nts_ax.set_xlabel('NTS', fontdict=LabelFont('label'))
	nts_ax.legend(loc='center left', bbox_to_anchor=(1.01,0.5), frameon=False, prop=LabelFont('legend'))
	ts_ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)
	nts_ax.grid(lw=1.2, ls='dashed', color='grey', zorder=1)

	return fig



def main(args=None):
	
	# loading params
	## preprocessing input files.
	samples = [transferPath(x) for x in args.input]
	sampleLabels = args.sampleLabels
	ref_label = args.ref_label
	scale_label = args.scale_label
	plotType = args.plotType
	outFig = transferPath(args.output)
	fig_title = args.title
	if not fig_title:
		fig_title = ''


	nSamples = len(samples)
	if not sampleLabels:
		sampleLabels = []
		for sm in samples:
			sampleLabels.append(os.path.basename(sm).split('.')[0])

	try:
		assert len(sampleLabels) == nSamples
	except AssertionError:
		info()
		sys.exit()


	# check if control group exists
	control_dict = None
	control = args.control
	if control:
		control = [transferPath(x) for x in control]
		control_dict = {}

		if len(control) == 1:
			_ctmp = _preprocessing_singlesample(control[0], ref_label, scale_label, plotType)
			for sm in samples:
				control_dict[sm] = copy.deepcopy(_ctmp)

			del _ctmp

		else:
			try:
				assert len(control) == nSamples
			except AssertionError:
				info()
				sys.exit()

			for idx, sm in enumerate(samples):
				control_dict[sm] = _preprocessing_singlesample(control[idx], ref_label, scale_label, plotType)


	# create dict to store infomation and data of samples
	# keys: sample_file_name, sample_mode, sample_mean_data, sample_tick_label, headerInfo, y_lable(method) ...
	samples_dict = {}
	for idx, sm in enumerate(samples):
		samples_dict[sm] = _preprocessing_singlesample(sm, ref_label, scale_label, plotType)
		samples_dict[sm]['sample_label'] = sampleLabels[idx]


	# preprocessing completed.
	# start to plot
	## set figure
	mergePlots = args.merge
	numPerRow = args.num_per_row
	#numRows = np.ceil(nSamples / numPerRow)
	numRows = int(np.ceil(nSamples / numPerRow))
	if nSamples <= numPerRow:
		numPerRow = nSamples

	if mergePlots:

		# Check all --input files are same in mode and length
		_mode_tmp = set()
		_cut_tmp = set()
		for k,v in samples_dict.items():
			sub_mode = v.get('mode')
			_mode_tmp.add(sub_mode)
			sub_cut = v.get('headerinfo')['SeparatorColumn']
			_cut_tmp.add(sub_cut)
		try:
			assert len(_mode_tmp) == 1 and len(_cut_tmp) == 1
		except AssertionError:
			info('choose merge but input file not in same in mode')
			sys.exit()

		#fig = plt.figure(constrained_layout=True)
		#gs = gridspec.GridSpec(1, 2, figure=fig)
		fig = merge_plot(samples_dict, control_dict)
		fig.suptitle(fig_title, fontdict=LabelFont('title'))
		fig.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()

	else:
		#fig = plt.figure(constrained_layout=True)
		#gs = gridspec.GridSpec(numRows, numPerRow, figure=fig)
		fig = split_plot(numRows, numPerRow, samples_dict, control_dict)
		fig.suptitle(fig_title, fontdict=LabelFont('title'))
		fig.savefig(outFig, dpi=300, bbox_inches='tight')
		plt.close()


if __name__ == '__main__':
	main()







###################
# To be modified
# 1: figsize
#    according to number of input files or subplots
# 2: suptitle font
##################