#! /usr/bin/env python2.6.6
"""
Reads an aligned SAM file and groups reads by position and then by molecular tag.
Input is aligned SAM file, output is aligned SAM file.
"""

## alphabetize -- stylistic but helpful
import argparse
import csv
import locale
import sets
import sys


def main(input_file, output_file):
    check_inputs(input_file, output_file)
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

    ## This "with" call is a common python idiom used with opening files and is called a "context manager."
    ## Context managers perform an action on opening and then another on closing, no matter what happens.
    ## When you use this "with" statment while opening a file, you can be sure that, no matter what,
    ## it will be closed when you exit. This is a good practice to adapt and prevents files from being left open
    ## in the event of an unexpected error.

    ## You usually only need the 'b' (e.g. "rb" instead of "r") for text files when you care about Windows:
    ## essentially it adds an additional end-of-line character as required by that OS.
    ## But anyway if you do it for one you may as well do it for both
    ## See http://stackoverflow.com/questions/15750660/python-file-io-w-vs-wb

    with open(input_file, 'rb') as in_handle, open(output_file, 'ab') as out_handle:
    ## more descriptive variable name helps keep track of things later ("sam_reader" vs. "reader")
        sam_reader = csv.reader(in_handle, delimiter = '\t')
        sam_writer = csv.writer(out_handle, delimiter = '\t')

        # Initialize variables
        ## Changed from "pos" to be a little more descriptive, just for others who read your code
        ## Also it's good to initialize "non-value" variables to None, which evaluates to False in python truth tests
        ## this is especially true if 0 or "" (the empty string) could be misinterpreted as a valid value
        pos_previous    = None
        ## another way to do this is "barcodes, lengths, quals, samrows = [], [], [], []
        ## but either is fine
        barcodes        = []
        lengths         = []
        quals           = []
        samrows         = []

        # Start reading sam file
        for row in sam_reader:

            ## one way to improve code clarity is to use descriptive variable names.
            ## a common way to do that is to expand lists into descriptive variables like this
            ## these names actually aren't great names but they're the ones from the SAM format document so whatever
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = row

            # Header lines start with @ and get written straight through to output file
            ## this "startswith" is a python string function, but qname[0] is fine too.
            if qname.startswith("@"):
                sam_writer(output_file, 'a', row)
                # no need to perform further checks on this line -- skip to the next "for" value
                next

            # If the read has the same starting position as the last one, append values
            ## You can move the first-time initialization to the "else" clause below because it will evaluate the same way as "new position" reads
            if pos == pos_previous:
                barcodes.append(qual)
                lengths.append(tlen)
                quals.append(mapq)
                samrows.append(row)
                # skip to next line
                next

            # If different starting position than last read, divide up according to length and reinitialize lists
            else:
                # Initialize for a new set of alignments
                pos_previous    = pos
                ## This is clearer and a bit cleaner-looking than like barcode = []; barcode.append(value), but either works
                barcodes        = [seq]
                lengths         = [tlen]
                quals           = [mapq]
                samrows         = [row]

            samrows = identify_subset(barcodes, lengths, qual, samrows)
            sam_writer.writerows(samrows)

    locale.setlocale(locale.LC_ALL, '')


## One thing you could do here is create a "sam_writer" class,
## which could be a fun exercise if you're interested:
## you could write the write() and writerows() functions yourself,
## And you'd need an __init__ of course too.
## otherwise I'd just use the csv.writer and call its functions (csv.writerows)
## You could actually create an object class that inherits from csv.writer even!
#def sam_write(output_file, mode, samrows):
#    """
#    Writes a list to a tab-delimited file
#    """
#    out_handle = open(output_file, mode)
#    writer = csv.writer(out_handle, delimiter = '\t')
#    for line in samrows:
#        writer.writerow(line)
#    out_handle.close()

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
