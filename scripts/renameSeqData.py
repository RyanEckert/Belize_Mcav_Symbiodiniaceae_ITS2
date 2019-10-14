#!/usr/bin/env python3

import os
import re

path = os.getcwd()
pattern1 = re.compile(r'S(\d)+_')
pattern2 = re.compile(r'RE-P(\d)-')
targets = []

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq')]:
	os.rename(fileName, fileName.replace('_001',''))

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq')]:
	os.rename(fileName, fileName.replace('L001_',''))

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq')]:	
	m1 = pattern1.search(fileName)
	if m1:
		os.rename(fileName, fileName.replace(m1.group(0),''))

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq')]:
	m2 = pattern2.match(fileName)
	if m2:
		os.rename(fileName, fileName.replace(m2.group(0),''))
print('DONE')
