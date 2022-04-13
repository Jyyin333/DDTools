#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from ddtools.utilities import *

def main(args=None):

	numbers = args.number

	sys.stdout.write('Max number is: {:.0f}'.format(
		max(numbers)))


if __name__ == '__main__':
	main()


