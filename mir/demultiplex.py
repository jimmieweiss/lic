#! /usr/bin/env python2.7

import gzip
import sys
from Bio import SeqIO
import HTSeq
import subprocess
from subprocess import Popen, PIPE

in_handle_one = open(sys.argv[1], "rU")
in_handle_two = open(sys.argv[2], "rU")
out_handle_one = open(sys.argv[3], "w")
out_handle_two = open(sys.argv[4], "w")
barcode = str(sys.argv[5])
#qual_cutoff = int(sys.argv[5])

part1 = ['grep', '-B 1', '-A 2', '^'+barcode, sys.argv[2]]
p = ['\\']
part2 = ['grep', '-v']
part3 = [sys.argv[4]]
print part1
print part2
print part3
p1 = Popen(part1)#, stdout=PIPE)
print 't'
output = p1.communicate()[0]
p2 = Popen(part2, stdin=p1.stdout, stdout=PIPE,)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p2.communicate()[0]
#part1 + ['>'] + part3
#print pipe
#output = subprocess.call(pipe)

#for record in HTSeq.FastqReader(in_handle_two):
#	if str(record.seq[:6]) == barcode:
#		record.write_to_fastq_file(out_handle_two)
in_handle_two.close()
out_handle_two.close()

out_handle_two = open(sys.argv[4], "rU")

parser1 = SeqIO.parse(in_handle_one, 'fastq')
#record2 = parser2.next()
for record2 in SeqIO.parse(out_handle_two, 'fastq'):
	for record1 in parser1:
		if record1.id == record2.id:
			#if sum(record1.letter_annotations["phred_quality"])/len(record1.seq) >= qual_cutoff:
			#	if sum(record2.letter_annotations["phred_quality"])/len(record2.seq) >= qual_cutoff:
			#record.write_to_fastq_file(out_handle_one)
			SeqIO.write(record1, out_handle_one, "fastq")
			break
		#record2 = parser2.next()

in_handle_one.close()
out_handle_one.close()
out_handle_two.close()


     
#def MoveBarcode(self):
		#print self.args.in_handle
#		for record in HTSeq.FastqReader(self.args.in_handle):
#			tag =  record.seq[:8]
#			record = record[8:]
#			record.name = record.name.split(' ')[0] + tag
#			record.write_to_fastq_file(self.args.out_handle)