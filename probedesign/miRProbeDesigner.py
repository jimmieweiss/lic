#! /usr/bin/env python2.7

import HTSeq
import sys
import pysam
import gzip
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.pyplot as plt
import pylab
import csv
import argparse

#
#   python miRProbeDesign.py probes MIRFILE 
#
#	skriver ut seed (forsta baserna med perfekt matchning) antal miR som fangas antal prober som behovs for att fanga dem

def main(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets, mode):
	
	if arg_vars['mode'] == True:
		print 'new mode'
		target_design(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets)
	
	elif arg_vars['perfect_matches']:
		probedesign(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets)
		#raise SyntaxError("Must indicate index length to tally.")
		#raise SyntaxError("Must indicate file to parse for indexes.")

def target_design(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets):
	fasta = FastaFile()
	fasta.read_in_mirlist()
	target_dict = fasta.return_targets()
	f = MiRList(target_dict)
	f.read_list()
	f.create_file()

def probedesign(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets):
	
	number_targets = number_targets - 2
	in_file = FastaFile(input_file, output_file, perfect_matches)
	in_file.Create_first_set()
	in_file.read_in_targets()
	probes = in_file.Create_first_probe_list(start_base)

	for key in probes:
		temp = in_file.Expand_probe(probes[key], mismatches, start_base, end_base)
		length = 0
		#print len(temp) #.keys()   
		while len(temp) > length:
			length = len(temp)
			temp = in_file.Expand_probe(temp, mismatches, start_base, end_base)
		targets, filtered = in_file.FilterOnTargets(temp, number_targets)
		#print filtered
		write_probes(output_file, filtered, targets, perfect_matches)
		#print temp   Ta bort # for att skriva alla ofiltrerade sekvenser och targets till stdout

def write_probes(output_file, filtered, targets, perfect_matches):

	if not arg_vars['output_file']:
		if len(targets) > 0:
			print str(filtered.keys()[0][0:perfect_matches]) + '  ' + str(len(targets)) + '  ' +str(len(filtered.keys()))
			for key in filtered.keys():
				print key + '    ' + str(filtered[key])   ### Ta bort # for att skriva alla sekvenser och deras target till stdout
		
	elif arg_vars['output_file']:
		if len(targets) > 0:
			print str(filtered.keys()[0][0:perfect_matches]) + '  ' + str(len(targets)) + '  ' +str(len(filtered.keys()))
			for key in filtered.keys():
				print key + '    ' + str(filtered[key])   ### Ta bort # for att skriva alla sekvenser och deras target till stdout
	

class FastaFile:
	def __init__(self, input_file, output_file, perfect_matches):
		self.in_handle = input_file
		self.out_handle = output_file
		self.bases = ['A', 'T', 'C', 'G']
		self.seed_length = perfect_matches ### Andra for att fa langre perfekt match i borjan
		self.seed_list = []
		self.targets = {}
		self.probes = {}
		
	def return_targets(self):
		return self.targets
	
	def read_in_mirlist(self):
		self.__in_handle = open(self.in_handle, "rU")
		for rec in SeqIO.parse(self.__in_handle, "fasta"):
			if rec.seq not in self.targets:
				self.targets[rec.id] = rec
		self.__in_handle.close()
		
	def read_in_targets(self):
		self.__in_handle = open(self.in_handle, "rU")
		for rec in SeqIO.parse(self.__in_handle, "fasta"):
			if rec.seq not in self.targets:
				self.targets[rec.seq.back_transcribe()] = rec.id
		self.__in_handle.close()
				

	def Create_first_set(self):
		if self.seed_length == 0:
			pass
		elif self.seed_length == 1:
			self.seed_list = self.bases
		else:
			self.seed_list = self.bases
			print range(self.seed_length-1)
			for i in range(self.seed_length-1):
				__templist = []
				for j in range(len(self.seed_list)):
					for m in range(len(self.bases)):
						__temp_bc = self.seed_list[j] + self.bases[m]
						__templist.append(__temp_bc)
				self.seed_list = __templist
	
	
	def Create_first_probe_list(self, start_base):
		for seed in self.seed_list:
			temp_dict = {}
			for rec in self.targets:
				if str(rec[start_base:(start_base+self.seed_length)]) == seed:
					if str(rec) not in temp_dict.keys():
						temp_dict[str(rec)[start_base:]] = [self.targets[rec]]
					else:
						temp_dict[str(rec)[start_base:]] = temp_dict[str(rec)] + [self.targets[rec]]
			self.probes[seed] = temp_dict
		return self.probes
	
	def Unique_probes(self):
		temp_dict = {}
		list = []
		for key in set.keys():
			if str(key) in list:
				list.append(str(key))
			else:
				list.append(str(kyey))
	
	def Expand_probe(self, probe_set, mismatches, start_base, end_base):
		iter = 0
		temp_dict = {}
		end = end_base - start_base ### Andra for att fa langre eller kortare region 16
		while iter < len(probe_set):
			keylist = list(probe_set.keys())
			startkey = keylist[iter]
			keylist = keylist[iter+1:len(probe_set)+1]
			for rec in keylist:
				temp_list = [startkey[0:end], rec[0:end]]
				temp_seq = ''
				for base in range(end):
					if startkey[base] == rec[base]:
						temp_seq = temp_seq + rec[base]
					else:
						temp_seq = temp_seq + 'N'
				temp_Seq = Seq(temp_seq, generic_dna)
				if temp_Seq.count('N') < (mismatches+1): ### Andra for att tillata fler eller farre N i regionen
					if temp_seq not in temp_dict.keys():
						temp_set = set(probe_set[startkey] + probe_set[rec])
						temp_dict[temp_seq] = list(temp_set)
					else:
						temp_set = set(temp_dict[temp_seq] + probe_set[startkey] + probe_set[rec])
						temp_dict[temp_seq] = list(temp_set)
				else:
					temp_list = [startkey, rec]
					for key in temp_list:
						if key not in temp_dict.keys():
							temp_dict[key] = probe_set[key]
						else:
							temp_set = set(temp_dict[key] + probe_set[key])
							temp_dict[key] = list(temp_set)
			iter = iter + 1
		#print temp_dict
		#print len(temp_dict)
		return temp_dict
			
	
	def FilterOnTargets(self, probeset, number_targets):
		targets = []
		filtered = {}
		endpos = number_targets  ###Andra for att krava fler miR per prob -1 = minst ett miR per prob, 0 = minst tva miR per prob osv.
		iter = 15
		while iter > endpos:
			for probe in probeset.keys():
				found = True
				if len(probeset[probe]) > iter:
					for target in probeset[probe]:
						if target not in targets:
							targets.append(target)
							found = False
					if found == False:
						filtered[probe] = probeset[probe]
			iter = iter - 1
		return targets, filtered
		#if len(targets) > 0:
		#	print str(filtered.keys()[0][0:self.seed_length]) + '  ' + str(len(targets)) + '  ' +str(len(filtered.keys()))
		#	for key in filtered.keys():
		#		print key + '    ' + str(filtered[key])   ### Ta bort # for att skriva alla sekvenser och deras target till stdout
		

	def FilterOnNrOfNs(self, probeset):
		targets = []
		filtered = {}

		for probe in probeset.keys():
			found = True
			if len(probeset[probe]) > iter:
				for target in probeset[probe]:
					if target not in targets:
						targets.append(target)
						found = False
				if found == False:
					filtered[probe] = probeset[probe]
			iter = iter - 1
		if len(targets) > 0:
			print str(filtered.keys()[0][0:self.seed_length]) + '  ' + str(len(targets)) + '  ' +str(len(filtered.keys()))
		print filtered   ### Ta bort # for att skriva alla sekvenser och deras target till stdout

class MiRList:
	
	def __init__(self, target_dict):
		self.in_handle = open(sys.argv[3], "rU")
		self.targets = target_dict
		self.mir_list = []
		self.out_handle = open(sys.argv[4], "w")
		
	def read_list(self):
		for line in self.in_handle:
			mir = line.strip('\n')
			if mir[-1] == 'p':
				mir = 'hsa-miR' + mir[7:]
				self.mir_list.append(mir)
			else:
				mir1 = 'hsa-miR' + mir[7:] + '-3p'
				mir2 = 'hsa-miR' + mir[7:] + '-5p'
				self.mir_list.append(mir1)
				self.mir_list.append(mir2)
	
	def create_file(self):
		for mir in self.mir_list:
			if mir in self.targets.keys():
				sequence = self.targets[mir]
				sequence.id = mir
				SeqIO.write(sequence, self.out_handle, 'fasta')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--input-file",
								help = "Alignment file in SAM format with molecular barcode as last column. Required")
	
	parser.add_argument("-o", "--output-file",
								help = "Output file for duplication marked reads. Required.")
	
	parser.add_argument("-s", "--start-base", type = int, 
								help = "Left most position for probe on target. Required.")
	
	parser.add_argument("-p", "--perfect-matches", type = int, 
								help = "Number of perfect matches at beginning of target. Required.")
	
	parser.add_argument("-m", "--mismatches", type = int, 
								help = "Number of perfect matches at beginning of target. Required.")

	parser.add_argument("-e", "--end-base", type = int, 
								help = "Number of perfect matches at beginning of target. Required.")

	parser.add_argument("-t", "--number-targets", type = int, 
								help = "Number of perfect matches at beginning of target. Required.")

	parser.add_argument("-l", "--mode", type = bool, 
								help = "Number of perfect matches at beginning of target. Required.")
     
	arg_vars = vars(parser.parse_args())
	
	locals().update(arg_vars)
	main(input_file, output_file, perfect_matches, start_base, mismatches, end_base, number_targets, mode)


