#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
from setuptools import setup, find_packages


def readme():
	fi = open('README.md')
	README = fi.read()
	fi.close()
	return README


setup(
	name='DDTools',
	version='1.0.0',
	author='Jiayong Yin',
	author_email='782515144@qq.com',
	description='Analysis tool suite specially for NGS Data based on End-Labeling strategy.',
	long_description=readme(),
	license='LICENSE.txt',
	packages=find_packages(),
	scripts=['bin/ddtools'],
	classifiers=[
		'License :: Free For Home Use',
		'License :: OSI Approved :: MIT License',
		'Intended Audience :: Science/Research',
		'Programming Language :: Python :: 3.7'
	],
	python_requires='>=3',
	install_requires=[
		"numpy >= 1.19.2",
		"pysam >= 0.14.0",
		"py2bit >= 0.3.0",
		"pyfaidx >= 0.5.9.1",
		"scipy >= 1.6.0",
		"matplotlib >= 3.3.0"
	],
	zip_safe=False
)


