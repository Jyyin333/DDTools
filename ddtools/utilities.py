#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import string, random
import time
import gzip



class SnakeBlock:

	def __init__(self, projectname, counter = -1, blocks=dict(), sortidx=dict(), collects=dict()):
		self.projectname = projectname
		self.counter = counter
		self._blocks = blocks
		self._sortidx = sortidx
		self._collects = collects

	def constructBlock(self, rulename, plist, collect=False):

		self.counter += 1

		rule_title = 'rule {}'.format(rulename)
		plist.insert(0, rule_title)

		flag = False
		res = ''''''

		for idx,k in enumerate(plist):

			if k == 'shell':
				flag = True

			if k in ['input','output','shell','params','threads','log']:
				k = '\t' + k + ':' + '\n'
				res += k
				continue

			if k.startswith('rule')  and not k.startswith('rules'):
				k = k + ':' + '\n'
				res += k
				continue

			else:

				if flag:
					k = '\t\t' + k + '\n'
					res += k
					continue
				else:
					k = '\t\t' + k + ',' + '\n'
					res += k
					continue


		self._sortidx[self.counter] = rulename
		self._blocks[rulename] = res
		self._collects[rulename] = collect


	def constructVar(self, vardict):
		var = ""

		fasta_line = 'GENOME_FA = "{}"'.format(vardict.get('fasta'))
		index_line = 'ALIGNER_INDEX = "{}"'.format(vardict.get('index'))
		picard_line = 'picard = "{}"'.format(vardict.get('picard'))
		mapq_line = 'MAPQ = {}'.format(vardict.get('mapq'))

		var = '\n'.join(map(str, [fasta_line, index_line, picard_line, mapq_line])) + '\n'

		self._varline = var


	def constructHeadAndTail(self):

		self._header = '''rule all:\n\tinput:\n\t\texpand("{sample}.Done", sample=SAMPLES)\n'''

		rulesKeep = []
		for k, v in self._collects.items():
			if v == True:
				rulesKeep.append(k)

		tailList = ['rules.{}.output,'.format(x) for x in rulesKeep]
		tailLine = 'rule Done:\n\tinput:\n'
		for r in tailList:
			tmp = '\t\t' + r + '\n'
			tailLine += tmp

		tailLine += '\toutput:\n\t\ttouch("{sample}.Done")\n'

		self._tail = tailLine


	def writeFile(self, snakefile):

		# write environment variables
		with open(snakefile, 'w') as fo:
			fo.write(self._varline)
			fo.write('\n\n')

			# write samplelist
			smpath = './initial/' + self.projectname + '.sampleList'
			try:
				assert os.path.exists(os.path.abspath(smpath))
			except:
				info("No sample list file in {}.".format(smpath))
				sys.exit()

			fo.write("SAMPLES = [line.strip() for line in open('{}')]\n\n".format(os.path.abspath(smpath)))

			# write configfile
			conf = './initial/' + self.projectname + '_project.json'
			fo.write("configfile:'{}'\n\n".format(os.path.abspath(conf)))

			# write header
			fo.write(self._header)
			fo.write('\n\n')

			# write steps according to index
			indices = sorted(self._sortidx.keys(), reverse=False)
			for i in range(len(indices)):
				idx = indices[i]
				rulename = self._sortidx.get(idx)
				fo.write(self._blocks.get(rulename))
				fo.write('\n\n')
			
			# write tail
			fo.write(self._tail)
			fo.write('\n\n')




def info(logging, tar='out'):

	called_time = localtime()
	logging = str(logging)

	if tar == 'out':
		sys.stdout.write('>>> {}\t{}\n'.format(called_time, logging))

	if tar == 'error':
		sys.stdout.write('Error:\n')
		sys.stdout.write('{}\t{}\n'.format(called_time, logging))

	if tar == 'warn':
		sys.stdout.write('... Warning:{}\n'.format(logging))
		


def FileIO(path, name, mode):

	filename = name
	tar = os.path.join(path, filename) if filename != None else transferPath(path)

	if filename != None and filename.startswith('temp'):

		while True:

			filefix = filename.split('_')[-1]
			rstring = get_random_string(8, 'let') + 'tmp'
			filename = filefix + '_' + rstring

			tar = os.path.join(path, filename)

			if not os.path.exists(tar):
				break


	if mode == 'dir':
		try:
			os.mkdir(tar)
		except FileExistsError as filerr:
			info('Fail to create directory {}, it has existed.'.format(tar), 'error')
			sys.exit()


	if mode == 'file':
		try:
			open(tar, 'r').close()
			info('Fail to create File {}, it has existed.'.format(tar), 'error')
			sys.exit()
		except FileNotFoundError:
			open(tar, 'w').close()
			os.remove(tar)

	if mode == 'check':
		pass

	return tar


def get_random_string(n, mode=None):

	if mode == 'let':
		res = ''.join(map(str,random.sample(string.ascii_letters, n)))

	if mode == 'dig':
		res = ''.join(map(str,random.sample(string.digits, n)))

	if mode == None:
		res = ''.join(map(str,random.sample(string.ascii_letters+string.digits,n)))


	return res



def localtime():
	t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
	return t
	# '2022-01-09 11:21:57' 



def transferPath(rpath):
	cmd = r'cd ~ && pwd'
	userpath = os.popen(cmd).read().strip()
	
	if rpath.startswith('~/'):
		tmp = rpath[2:]
		res = os.path.join(userpath, tmp)
		#res = rpath.replace('~', userpath)
		return res
	
	if rpath.startswith('/'):
		return rpath

	if rpath.startswith('.'):
		res = os.path.abspath(rpath)
		return res
	
	else:
		pwd = os.popen('pwd').read().strip()
		return os.path.join(pwd, rpath)

def checkUnixTools(tool):
	
	toolname = tool
	userEnvirPath = os.popen('echo $PATH').read().strip().split(':')[:-1]
	flag = False
	for p in userEnvirPath:
		sp = os.path.join(p, toolname)
		if os.path.exists(sp):
			flag = True
			break

	return flag


def fetch_sequence(fasta, chrom, start, end):
	try:
		sequence = fasta[chrom][start:end]
	except KeyError:
		info("warning: {name} not found in file\n".format(**locals()), 'error')
		return
	except ValueError as ve:
		sys.stderr.write(str(ve))
		return
	return sequence


def writeHeader(target_file, header_dict):

	import json

	if target_file.endswith('.gz') or target_file.endswith('.bgz'):
		# step 1: read raw content
		try:
			with gzip.open(target_file, 'rb') as fi:
				raw_content = fi.read()

			# step 2: create a new file to write
			_tmpfilename = os.path.join(os.path.dirname(target_file), get_random_string(15)+'.gz')
			_tmpfile = gzip.open(_tmpfilename, 'wb')
			_tmpfile.write('#'.encode('utf-8'))
			_tmpfile.write(json.dumps(header_dict).encode('utf-8'))
			_tmpfile.write('\n'.encode('utf-8'))
			_tmpfile.write(raw_content)
			_tmpfile.close()

			# step 3: rename
			os.remove(target_file)
			os.rename(_tmpfilename, target_file)

		except FileNotFoundError:
			with gzip.open(target_file, 'wb') as fo:
				fo.write('#'.encode('utf-8'))
				fo.write(json.dumps(header_dict).encode('utf-8'))
				fo.write('\n'.encode('utf-8'))

	else:
		try:
			with open(target_file, 'r+') as fo:
				raw_content = fo.read()
				fo.seek(0)
				fo.write('#')
				fo.write(json.dumps(header_dict))
				fo.write('\n')
				fo.write(raw_content)

		except FileNotFoundError:
			with open(target_file, 'w') as fo:
				fo.write('#')
				fo.write(json.dumps(header_dict))
				fo.write('\n')


def readHeader(target_file):

	import json

	if target_file.endswith('.gz') or target_file.endswith('.bgz'):
		with gzip.open(target_file, 'rb') as fi:
			header = fi.readline().decode()

	else:
		with open(target_file, 'r') as fi:
			header = fi.readline()


	loads = header[1:]
	res = json.loads(loads)

	return res
	

def NormalizeParser():
	import argparse

	NormalizingParser = argparse.ArgumentParser('normalize',add_help=False)
	NormalizingParser.add_argument('--method', default='RPKM', choices=['RPKM','CPM','DEPTH'], help='Method to normalizing signal. (default: %(default)s)')
	NormalizingParser.add_argument('--readLength', type=int, help='The sequencing read length in bp.')
	NormalizingParser.add_argument('--genomeSize', help='Effective genome sizes')

	return NormalizingParser


def get_base_content(chrom, start, end, tbit, fraction=True, base = 'G'):

	# check if base valid
	base = base.upper().strip()
	for i in range(len(base)):
		try:
			assert base[i] in ['A','C','G','T']
		except AssertionError:
			sys.stderr.write('Invalid base type.\n')
			return None

	try:
		bases = tbit.bases(chrom, start, end, fraction=False)
	except RuntimeError:
		return 0
	except:
		return 0

	if end > tbit.chroms(chrom):
		end = tbit.chroms(chrom)
	if sum(bases.values()) < 0.95 * (start - end):
		raise Exception("WARNING: too many NNNs present in {}:{}-{}".format(chrom, start, end))
		return None

	if len(base) == 1:
		if fraction:
			return (bases[base]) / float(end - start)
		return bases[base]

	else:
		if fraction:
			return sum([bases[base[i]] for i in range(len(base))]) / float(end - start)
		return sum([bases[base[i]] for i in range(len(base))])


def RGB2hexadecimal(rgb):

	res = '#'
	rgb = rgb.split(',')
	for c in rgb:
		c = int(c)
		res += str(hex(c))[-2:].replace('x','0').upper()

	return res
	
def LabelFont(query):

	store = {
		'title':{
			'family':'Arial',
			'weight':'700',
			'size':16
		},
		'label':{
			'family':'Arial',
			'weight':'600',
			'size':15
		},
		'tick':{
			'family':'Arial',
			'weight':'500',
			'size':13
		},
		'legend':{
			'family':'Arial',
			'weight':'600',
			'size':14.5
		}
	}

	return(store.get(query))


def DefaultColors(idx):

	loop_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
	'#653700', '#ff81c0', '#bf77f6', '#ae7181', '#033500',
	'#6e750e', '#ad8150', '#06c2ac', '#cea2fd', '#cf6275',
	'#0165fc', '#ceb301', '#9aae07', '#fdaa48', '#ff474c']

	return loop_colors[idx]


def bedCountsinRegion(tbx, chrom, start, end, totalCounts=None, normlizeMethod=None, readLength=None, genomeSize=None):

	p, m = 0, 0
	tile_length = abs(end - start)
	try:
		records = tbx.fetch(chrom, start, end)
	except ValueError as ve:
		info(ve,'error')
		return 0, 0

	for read in records:

		if read.strand == '+':
			p += 1

		else:
			m += 1

	if normlizeMethod == 'RPKM':

		p = round(p / totalCounts / tile_length * 1e9, 6)
		m = round(m / totalCounts / tile_length * 1e9, 6)

	elif normlizeMethod == 'CPM':

		p = round(p / totalCounts * 1e6, 6)
		m = round(m / totalCounts * 1e6, 6)

	elif normlizeMethod == 'DEPTH':

		scale_factor = totalCounts * readLength / genomeSize
		p /= scale_factor
		m /= scale_factor

	# warnig: one possible return is raw read counts.

	return p, m


def process_refsite(frag_start,frag_end,frag_strand,rp='TSS'):

	ref_site = 0

	if rp == 'TSS':
		if frag_strand == '+':
			ref_site = frag_start
		if frag_strand == '-':
			ref_site = frag_end
		else:
			ref_site = frag_start

	elif rp == 'TES':
		if frag_strand == '+':
			ref_site = frag_end
		if frag_strand == '-':
			ref_site = frag_start
		else:
			ref_site = frag_end

	elif rp == 'CENTER':
		ref_site = int(round((frag_start+frag_end)/2,0))

	return ref_site


def process_regions(params_dict, mode):

	params_dict = params_dict.copy()
	regionfile = params_dict.get('regionfile')

	if regionfile.endswith('.gz') or regionfile.endswith('.bgz'):
		region_records = [line.decode().strip() for line in gzip.open(regionfile)]
	else:
		region_records = [line.strip() for line in open(regionfile)]


	upstream = params_dict.get('upstream')
	downstream = params_dict.get('downstream')

	res = []
	totalCounts = 0
	nFaileds = 0

	if mode == 'refpoint':

		rp = params_dict.get('rp')
		#binsize = params_dict.get('binsize')
		#stepsize = params_dict.get('stepsize')

		for line in region_records:

			totalCounts += 1

			if line.startswith('#') or len(line) == 0:
				continue

			line = line.split('\t')

			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			try:
				strand = line[5]
			except IndexError:
				strand = '.'

			ref_site = process_refsite(start,end,strand,rp=rp)

			if strand == '+' and ref_site - upstream < 0:
				nFaileds += 1
				continue

			if strand == '-' and ref_site - downstream < 0:
				nFaileds += 1
				continue

			if strand == '.' and ref_site - upstream < 0:
				continue

			res.append([chrom, ref_site, strand])


	if mode == 'scale':
		
		for line in region_records:

			totalCounts += 1

			if line.startswith('#') or len(line) == 0:
				continue

			line = line.split('\t')

			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			try:
				strand = line[5]
			except IndexError:
				strand = '.'

			if strand == '+':
				up_region = [start - upstream, start]
				mid_region = [start, end]
				down_region = [end, end + downstream]

				if up_region[0] < 0:  # down_region[1] > chrom_max_size; but we cannot acquire the infomation about the maximum chromosome length in here. Therefore this is a potential bug.
					nFaileds += 1
					continue

			if strand == '-':
				up_region = [end, end + upstream]
				mid_region = [start, end]
				down_region = [start - downstream, start]

				if down_region[0] < 0:
					nFaileds += 1
					continue

			res.append([chrom, strand, up_region, mid_region, down_region])


	filter_info = 'Total regions: {:.0f}'.format(totalCounts)
	if nFaileds > 0:
		filter_info = 'Of total {:.0f} input regions, {:.0f} were discarded because of out of range.'.format(totalCounts, nFaileds)

	return res, filter_info


def filterBed(bedfile, ref_fa, output, keep_bases):

	base_length = list(set([len(x) for x in keep_bases]))

	assert len(base_length) == 1

	fo = open(output, 'w')

	if base_length[0] == 1:

		if bedfile.endswith('.gz') or bedfile.endswith('.bgz'):

			with gzip.open(bedfile, 'r') as fi:

				for line in fi:

					line = line.decode().strip().split('\t')

					if line[3] not in keep_bases:
						continue

					else:
						ol = '\t'.join(map(str, line)) + '\n'
						fo.write(ol)

		else:
			with open(bedfile, 'r') as fi:

				for line in fi:

					line = line.strip().split('\t')

					if line[3] not in keep_bases:
						continue

					else:
						ol = '\t'.join(map(str, line)) + '\n'
						fo.write(ol)

	elif base_length[0] == 2:

		#ref_fa = Fasta(transferPath(fasta))

		if bedfile.endswith('.gz') or bedfile.endswith('.bgz'):

			with gzip.open(bedfile, 'r') as fi:

				for line in fi:

					line = line.decode().strip().split('\t')

					chrom = line[0]
					start = int(line[1])
					end  = int(line[2])
					strand = line[5]

					if strand == '+':
						di_base = fetch_sequence(ref_fa, chrom, start-1, end)
						di_base = di_base.seq.upper()
					else:
						di_base = fetch_sequence(ref_fa, chrom, start, end+1)
						di_base = di_base.reverse.complement.seq.upper()

					if di_base not in keep_bases:
						continue
					else:
						ol = '\t'.join(map(str, line)) + '\n'
						fo.write(ol)

		else:

			with open(bedfile, 'r') as fi:

				for line in fi:

					line = line.strip().split('\t')

					chrom = line[0]
					start = int(line[1])
					end  = int(line[2])
					strand = line[5]

					if strand == '+':
						di_base = fetch_sequence(ref_fa, chrom, start-1, end)
						di_base = di_base.seq.upper()
					else:
						di_base = fetch_sequence(ref_fa, chrom, start, end+1)
						di_base = di_base.reverse.complement.seq.upper()

					if di_base not in keep_bases:
						continue
					else:
						ol = '\t'.join(map(str, line)) + '\n'
						fo.write(ol)


	fo.close()
	ref_fa.close()

	# creat bed index
	bgzip_cmd = 'bgzip {}'.format(output)
	bgziped_file = output + '.gz'
	tabix_cmd = 'tabix -p bed {}'.format(bgziped_file)

	os.system(bgzip_cmd)
	#time.sleep(15) # need to optimize

	########
	# optimize
	latest_size = 1
	zipfile_size = -1
	while latest_size > zipfile_size:
		zipfile_size = latest_size
		latest_size = os.path.getsize(bgziped_file)
		time.sleep(2)

	time.sleep(5)
	os.system(tabix_cmd)
	time.sleep(10)

	tmpfile_list = [bgziped_file, output+'.gz.tbi']

	return tmpfile_list



def checkIndex(path):
	
	# hg19.amb  hg19.ann  hg19.bwt  hg19.pac  hg19.sa
	# hg19.1.bt2  hg19.2.bt2  hg19.3.bt2  hg19.4.bt2  hg19.rev.1.bt2  hg19.rev.2.bt2
	bwa_suffixs = ['.amb', '.ann', '.bwt', '.pac', '.sa']
	bowtie2_suffixs = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

	bwa_search = set([os.path.exists(path+x) for x in bwa_suffixs])
	bt2_search = set([os.path.exists(path+x) for x in bowtie2_suffixs])

	if True in bwa_search | bt2_search:
		return True

	return False




		
				