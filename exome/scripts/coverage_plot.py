#! /usr/bin/env python2.7

import gzip
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import HTSeq
import math

in_handle = open('coveredbykit/coveredbykit.sample_interval_summary')
reader = csv.reader(in_handle, delimiter = '\t')
depth= [0]
s=0
chr='chr1'
for line in reader:
	if s == 0:
		s =1
	elif s < 1000:
		#pos.append(line[1])
		depth.append(float(line[2]))
		#print m
		s = s+1

length = len(depth)#[-1]+1
#print length

l = np.zeros(length)

maxx= max(depth)
for i in range(len(depth)):
	#print l[i]
	depth[i] = float(depth[i])/2+100#/(maxx) +1
tt = np.arange(0, 2*np.pi, 2*np.pi/len(depth))

d = []
for i in range(len(depth)):
	#print l[i]
	d.append(float(depth[-i])/2+180)

ratio = []
for i in range(len(depth)):
	#print l[i]
	ratio.append(((depth[i]-100)-(d[i]-300))/4+300)
#print tt
#print depth

ax = plt.subplot(111, polar=True)
#ax.plot(cr1, d, color='g', linewidth=3)
ax.plot(tt, depth, color='g', linewidth=1)
ax.plot(tt, d, color='b', linewidth=1)
ax.plot(tt, ratio, color='r', linewidth=1)
plt.axis('off')
#ax.set_rmax(0.0001)
#ax.set_rmax(max(d)+1)
ax.set_title('Coverage')
plt.show()
in_handle.close()