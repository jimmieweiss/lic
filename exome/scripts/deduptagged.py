#! /usr/bin/env python2.6.6
"""
Reads an aligned SAM file and groups reads by position and then by molecular tag.
Input is aligned SAM file, output is aligned SAM file.
"""

from __future__ import print_function

## I usually organize the imports alphabetically
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
    ## This is a new-ish print function. To use it, you need to do the "from __future__" bit at the top
    ## I prefer it and it's recommended because there's no alternative in python 3, which is the future of python.
    ## See that you can specify the output file for the print statement, and I've chosen STDERR (on the command line 1 = STDOUT, 2 = STDERR).
    ## In this case, sys.stdout would be fine, but often programs will print useful output to STDOUT and error to STDERR, so that the user
    ## can redirect STDOUT to a file but still see the errors printed to the screen (or redirect them elsewhere). For example,
    ## a user could do python deduptagged.py -i input_file.txt 1> output.txt 2>output.err
    print("Reading data from file \"{}\"...".format(input_file), end="", file=sys.stderr)
    total_data_lines        = 0
    ## python's built-in collections module is great. The deque object is like a list except it works fastest
    ## when you're dealing with the edges of it -- like popping things off the left or right side.
    ## Technically I think it's implemented as a "linked list," which you can read about if you like.
    ## My favorite collections object is defaultdict, though. It's awesome. You should read about that even if you don't like.
    non_duplicate_indexes   = collections.deque()

    ## My approach wound being to basically to run through the file and note the highest-quality read for each subgroup.
    ## This took me a little while to think of and is essentially what you were doing; however, the biggest difference
    ## is that it takes an iterative approach -- it doesn't deal with the entire list all at once,
    ## it deals with one element at a time. Python is built so that iterative approaches are easy, and in this case
    ## (and cases like it where you're dealing with really big files) it can be much faster because sorting and finding
    ## elements in really big lists takes a very very long time.
    ## I think if you approach this problem with a speed-related mindset, this iterative approach comes sort of naturally,
    ## but it's kind of weird to think this way. Anyway, when you're parsing 40GB files, every function call counts.

    ## Another way to think of it is, "What do I want to do here? I want to sort the reads by position, length, and barcode,
    ## and then find the member of that group with the higest quality, and mark all the others as duplicates."

    ## This "with-as" call is extremely Pythonic. It takes advantage of a Python feature called "context managers."
    ## Context managers perform an action when they start and another action when they're finished;
    ## in this case, it opens the file when it starts and it closes it when it's finished.
    ## This is better in general than manually opening and closing a file because a) it's faster to open and close once,
    ## and b) if there's some error in the middle it still closes the file (whereas a close() statement might never be reached).
    with open(input_file, 'rb') as in_handle:
        sam_reader = csv.reader(in_handle, delimiter = '\t')
        max_quals_dict = {}
        ## Enumerate is a function that works on all iterators (lists, tuples, strings even).
        ## It returns a tuple of (index, value) and is more Pythonic (and easier!!)
        ## than the corresponding "for i in range(len([0,1,2,3]))"
        for index, row in enumerate(sam_reader):
            if row[0][0] == ("@"):
                # Skip header lines
                next
            else:
                ## Using descriptive variable names is better than indexes for other people reading your code
                ## (which includes you sometimes, months after you wrote it)
                position    = int(row[3])
                qual        = int(row[4])
                length      = int(row[8])
                barcode     = row[-1]
                total_data_lines += 1
                ## This is kind of a complicated if statement, but basically it reads:
                ## If there is not yet a quality value for this combination of position, barcode, and length,
                ## or if it does exist and the current quality is higher than that one,
                ## update it with this value (and note the position of this read in the file).
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
                    ## This is a little weird too I guess. Basically it checks the list to see if the current index
                    ## was "marked" as a non-duplicate, and if so, it pops that off the list and
                    ## continues without marking the actual read data.
                    non_duplicate_indexes.remove(index)
                    # This is not a duplicate
                except ValueError:
                    ## However, if this one is -not- marked as a non-duplicate, it adjusts the appropriate flag before writing.
                    # This is a duplicate
                    row[1] += str(int(row[1]) + 0x400)
                ## Everything gets written -- the difference is just if their flags were modified or not
                sam_writer.writerow(row)
    print(" complete.", file=sys.stderr)
    duplicate_reads     = total_data_lines - non_duplicate_reads
    duplication_rate    = float(duplicate_reads) / total_data_lines
    ## String formatting is great
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
    # Note that csv.reader and csv.writer can take any iterator (files are iterators too!) as the source/destination of data!
    # This includes filehandles, lists, and even things like sys.stdout, sys.stdin
    #parser.add_argument("-s", "--use-stdbuf",
    #                            help="Read data from STDIN, write to STDOUT.")

    arg_vars = vars(parser.parse_args())
    locals().update(arg_vars)
    main(input_file, output_file)
