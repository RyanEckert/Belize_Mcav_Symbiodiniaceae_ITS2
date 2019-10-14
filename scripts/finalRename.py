#!/usr/bin/env python3

import os
import re

path = os.getcwd()
pattern = re.compile(r'.final')
targets = []

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq')]:	
	m = pattern.search(fileName)
	if m:
		os.rename(fileName, fileName.replace(m.group(0),''))

print('DONE')
