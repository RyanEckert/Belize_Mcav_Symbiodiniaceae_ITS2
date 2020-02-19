#!/usr/bin/env python3

import os, csv, argparse, sys


parser = argparse.ArgumentParser(description='Rename files based on columns of data from a specified .csv file')

parser.add_argument("-f", type=str, help="file extension type (not including '.'). [default = fastq]", action = "store", default='fastq')

parser.add_argument("-i", type=str, help="input .csv file name (not including '.csv') Should include filenames and sample numbers as columns, and be located in the working directory. [default = sampleList]", action = "store", default='sampleList')

parser.add_argument("-n", type=str, help="common name given to all samples, preceding the sample number in .csv file, if used.", action = "store", default='')

args = parser.parse_args()
i = args.i
f = args.f
n = args.n

with open(i + '.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')

    for row in csvreader:
        name = row[0]+'.'+f
        new = n+row[1]+'.'+f
        if os.path.exists(name):
            os.rename(name, new)

        else:
            print(name + " does not exist")
