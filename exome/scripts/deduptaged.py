#! /usr/bin/env python2.6.6

import sys
import csv
import sets
import locale
import argparse


def main(input_file, output_file):
	check_inputs(input_file, output_file)
	locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
	in_handle = open(input_file, 'rb')
	reader = csv.reader(in_handle, delimiter = '\t')
	
	
	# Bunch of variables
	pos = 0
	barcodes = []
	lengths = []
	qual = []
	samrows = []
	
	# Start reading sam file
	for row in reader:
		
		# Read header and write to output file
		if row[0][0] == '@':
			out_handle = open(output_file, 'a')
			writer = csv.writer(out_handle, delimiter = '\t')
			writer.writerow(row)
			out_handle.close()
			#sam_writer(output_file, 'a', row)
			
		# If same starting position as last read
		elif row[3] == pos:
			barcodes.append(row[-1])
			lengths.append(row[8])
			qual.append(row[4])
			samrows.append(row)
			
	    # First iteration
		elif pos == 0:
			barcodes.append(row[-1])
			lengths.append(row[8])
			pos = row[3]
			qual.append(row[4])
			samrows.append(row)
			
	    # If different starting position than last read
    	# Start dividing up according to length
		else:
			samrows = identify_subset(barcodes, lengths, qual, samrows)
			
			# Write alignments to output file
			sam_writer(output_file, 'a', samrows)
				
			# Clear history for new set of alignments
			pos = row[3]
			barcodes = []
			barcodes.append(row[-1])
			lengths = []
			lengths.append(row[8])
			qual = []
			qual.append(row[4])
			samrows = []
			samrows.append(row)
	
	# Take care of last set of lines	
	samrows = identify_subset(barcodes, lengths, qual, samrows)
	sam_writer(output_file, 'a', samrows)
	
	in_handle.close()
	locale.setlocale(locale.LC_ALL, '')
	
	
# SAM writer	
def sam_writer(output_file, mode, samrows):	
	out_handle = open(output_file, mode)
	writer = csv.writer(out_handle, delimiter = '\t')
	for line in samrows:
		writer.writerow(line)
	out_handle.close()

# Check validity of input variables.
def check_inputs(input_file, output_file):
	if not input_file:
		raise SyntaxError("Must indicate an input file.")
	if not output_file:
		raise SyntaxError("Must indicate an output file.")
	else:
		pass

# Divide set according to length and process them as a subset
def identify_subset(barcodes, lengths, qual, samrows):
	for length in set(lengths):
		subbar = []
		subqual = []
		subsam = []
		
		# Collect data for subset of barcodes
		for i in range(len(barcodes)):
			if lengths[i] == length:
				#if int(qual[i])>= 1:
				subbar.append(barcodes[i])
				subqual.append(qual[i])
				subsam.append(samrows[i])
				
		# If subset of reads > 0 -> find duplicates
		if len(subbar) > 0:
			samrows = identify_duplicates(subbar, subsam, subqual, samrows)
	return samrows

# Identify duplicates and mark all but the one with the best alignment.
def identify_duplicates(subbar, subsam, subqual, samrows):

	# Check each molecular barcode for duplicates
	for item in set(subbar):
		if subbar.count(item) > 1:
			dupcand = []
			dupcandqual = []
			
			# Collect data for duplicate candidates
			for j in range(len(subbar)):
				if subbar[j] == item:
					dupcand.append(subsam[j])
					dupcandqual.append(int(subqual[j]))
					
			# Choose one read not to marked
			notdup = dupcand[dupcandqual.index(max(dupcandqual))]
			
			# Add duplication to the SAM flag
			for dup in dupcand:
				if dup != notdup:
					samdup = samrows.index(dup)
					samrows[samdup][1] = int(samrows[samdup][1]) + 1024
	return samrows

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--input-file",
								help="Alignment file in SAM format with molecular barcode as last column. Required")
	
	parser.add_argument("-o", "--output-file",
								help="Output file for duplication marked reads. Required.")
    
	arg_vars = vars(parser.parse_args())
	locals().update(arg_vars)
	main(input_file, output_file)
