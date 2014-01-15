#! /usr/bin/env python2.6.6
"""
Reads an aligned SAM file and groups reads by position and then by molecular tag.
Input is aligned SAM file, output is aligned SAM file.
"""

from __future__ import print_function

import argparse
import collections
import csv
import locale
import sets
import sys

# TODO implement reading, writing to stdin,stdout
# TODO implement reading into / writing from memory, instead of reading the file twice
#      --> requires checking available memory, or trusting user to pass the appropriate flag
def main(input_file, output_file):
    check_inputs(input_file, output_file)
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    print("Reading data from file \"{}\"...".format(input_file), end="", file=sys.stderr)
    total_data_lines        = 0
    non_duplicate_indexes   = collections.deque()
    with open(input_file, 'rb') as in_handle:
        sam_reader = csv.reader(in_handle, delimiter = '\t')
        max_quals_dict          = {}
        # Start reading sam file
        for index, row in enumerate(sam_reader):
            # Header lines start with @ and get written straight through to output file
            if row[0][0] == ("@"):
                # Skip header lines
                next
            else:
                position    = int(row[3])
                qual        = int(row[4])
                length      = int(row[8])
                barcode     = row[-1]
                total_data_lines += 1
                if (not max_quals_dict.get((position, barcode, length))) \
                    or qual > max_quals_dict.get((position, barcode, length)):
                    # This is the highest-quality element at this position of this length with this molecular barcode thus far
                    max_quals_dict[(position, barcode, length)] = qual
                    non_duplicate_indexes.append(index)
    print(" complete.", file=sys.stderr)
    print("Writing data to file \"{}\"...".format(output_file), end="", file=sys.stderr)
    # The other way to do this is to load everything into memory -- probably faster but wouldn't work for huge files(?)
    non_duplicate_reads = len(non_duplicate_indexes)
    with open(input_file, 'rb') as in_handle, open(output_file, 'wb') as out_handle:
        sam_reader = csv.reader(in_handle, delimiter = '\t')
        sam_writer = csv.writer(out_handle, delimiter= '\t')
        for index, row in enumerate(sam_reader):
            if row[0][0] == ("@"):
                # Write headers straight through
                sam_writer.writerow(row)
            else:
                try:
                    # This may not be the fastest way to do this
                    non_duplicate_indexes.remove(index)
                    # This is not a duplicate
                except ValueError:
                    # This is a duplicate
                    # 1024 or 1023? I'm unclear on this
                    row[1] += str(int(row[1]) + 0x400)
                finally:
                    sam_writer.writerow(row)
    print(" complete.", file=sys.stderr)
    duplicate_reads     = total_data_lines - non_duplicate_reads
    duplication_rate    = float(duplicate_reads) / total_data_lines
    print(  "Wrote {total_data_lines} data lines:\n\t" \
            "{duplication_rate:0.2f}% duplication rate ({duplicate_reads} reads)".format(
                total_data_lines    = total_data_lines,
                duplication_rate    = 100 * duplication_rate,
                duplicate_reads     = duplicate_reads,), file=sys.stdout)
    locale.setlocale(locale.LC_ALL, '')


def check_inputs(input_file, output_file):
    """
    Check validity of input variables.
    """
    if not input_file:
        raise SyntaxError("Must indicate an input file.")
    if not output_file:
        raise SyntaxError("Must indicate an output file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input-file",
                                help="Alignment file in SAM format with molecular barcode as last column. Required")

    parser.add_argument("-o", "--output-file",
                                help="Output file for duplication marked reads. Required.")
    # TODO implement this!
    #parser.add_argument("-s", "--use-stdbuf",
    #                            help="Read data from STDIN, write to STDOUT.")

    arg_vars = vars(parser.parse_args())
    locals().update(arg_vars)
    main(input_file, output_file)
