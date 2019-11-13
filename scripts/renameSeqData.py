#!/usr/bin/env python3

import os
import re

path = os.getcwd()
pattern1 = re.compile(r'S(\d)+_L(\d)+_')
targets = []

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq.gz')]:  
        m1 = pattern1.search(fileName)
        if m1:
                os.rename(fileName, fileName.replace(m1.group(0),''))

for fileName in [fileName for fileName in os.listdir(path) if fileName.endswith('.fastq.gz')]:
        os.rename(fileName, fileName.replace('_001.','.'))

print('DONE')
