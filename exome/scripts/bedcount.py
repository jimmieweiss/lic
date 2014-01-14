#! /usr/bin/env python2.7

import gzip
import sys
import csv

reader = csv.reader(open(sys.argv[1], "rb"), delimiter='\t')
writer = csv.writer(open(sys.argv[2], "w"), delimiter='\t') 
basesum = 0
for line in reader:
	try:
		bases = int(line[2]) - int(line[1])
		basesum = basesum + bases
		line.append(bases)
	except:
		pass
	writer.writerow(line)

print basesum