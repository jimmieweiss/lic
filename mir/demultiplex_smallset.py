#! /usr/bin/env python2.7

import gzip
import sys
from Bio import SeqIO
#import HTSeq
import subprocess
#from subprocess import Popen, PIPE

in_handle_one = open(sys.argv[1], "rU")
in_handle_two = open(sys.argv[2], "rU")
out_handle = open(sys.argv[3], "w")

print 'start'

id_list=[]
for record in SeqIO.parse(in_handle_one, 'fastq'):
	id_list.append(record.id)
print 'id list done'
print id_list[0]
counter = 0
for record in SeqIO.parse(in_handle_two, 'fastq'):
		if record.id == id_list[counter]:
			SeqIO.write(record, out_handle, "fastq")
			print counter
			counter = counter + 1
			

in_handle_one.close()
in_handle_two.close()
out_handle.close()

