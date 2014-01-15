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

    #with open(input_file, 'rb') as in_handle, open(output_file, 'ab') as out_handle:
    with open(input_file, 'rb') as in_handle:
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
        barcodes_list   = []
        lengths_list    = []
        quals_list      = []
        samrows_list    = []

        #max_quals_dict = collections.defaultdict(dict)
        max_quals_dict          = {}
        non_duplicate_indexes   = []

        # Start reading sam file
        for index, row in enumerate(sam_reader):

            # Header lines start with @ and get written straight through to output file
            if row[0][0] == ("@"):
                # Skip header lines
                next
            else:
                ## one way to improve code clarity is to use descriptive variable names.
                position    = row[3]
                qual        = row[4]
                length      = row[8]
                barcode     = row[-1]

                if max_quals_dict.get((position, barcode, length)) and qual > high_quals_dict.get((position, barcode, length)).get(qual):
                    # This is the highest-quality element at this position of this length with this molecular barcode thus far
                    max_quals_dict((position, barcode, length)) = qual
                    #max_quals_dict((position, barcode, length))["index"]    = index
                    non_duplicate_indexes.append(index)

    # The other way to do this is to load everything into memory -- probably faster but wouldn't work for huge files(?)
    with open(input_file, 'rb') as in_handle, open(output_file, 'ab') as out_handle:
        for index, row in enumerate(sam_reader):
            if row[0][0] == ("@"):
                # Write headers straight through
                sam_writer.writerow(row)
            else:
                try:
                    l1.remove(index)
                    # This is not a duplicate
                except ValueError:
                    # This is a duplicate
                    # 1024 or 1023? I'm unclear on this
                    row[1] += 0x400
                finally:
                    samwriter.writerow(row)

    sys.exit()

#                # If the read has the same starting position as the last one, append values
#                ## You can move the first-time initialization to the "else" clause below because it will evaluate the same way as "new position" reads
#                if pos == pos_previous:
#                    barcodes_list.append(barcode)
#                    lengths_list.append(length)
#                    quals_list.append(qual)
#                    samrows_list.append(row)
#                    next
#                # If different starting position than last read, divide up according to length and reinitialize lists
#                else:
#                    # TODO change this it's ugly
#                    if pos_previous:
#                        # Interestingly, when you pass a list in python you pass the reference to the real object
#                        # This means that anything you do to the list within that function happens to the actual list
#                        # (This is true for all "mutable types" -->  see http://www.spontaneoussymmetry.com/blog/archives/438)
#                        # This is good to know as it's different than e.g. perl or C++, which pass copies of objects (I think?)
#                        # Anyway long story short, you don't need to reassign samrows because you can edit it within the function
#                        mark_duplicates(barcodes_list, lengths_list, quals_list, samrows_list)
#                        #samrows = identify_subset(barcodes, lengths, qual, samrows)
#                        sam_writer.writerows(samrows)
#                    # Initialize for a new set of alignments
#                    pos_previous    = pos
#                    ## This is clearer and a bit cleaner-looking than like barcode = []; barcode.append(value), but either works
#                    barcodes        = [barcode]
#                    lengths         = [length]
#                    quals           = [qual]
#                    samrows         = [row]
#
        #mark_duplicates(barcodes_list, lengths_list, quals_list, samrows_list)
        #samrows = identify_subset(barcodes, lengths, qual, samrows)
        #sam_writer.writerows(samrows)

    #locale.setlocale(locale.LC_ALL, '')


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
def identify_subset(barcodes, lengths, quals, samrows):
    import pdb; pdb.set_trace()
    ## Instead of looping through the list multiple times (one time per length)
    ## Use a dict (hash) where the keys are the length
    ## And the values are lists of all the matching rows
    #for length in set(lengths):
    rows_by_length = collections.defaultdict(list)
    ## collections are awesome datatypes for python -- defaultdict is a good one: http://blog.ludovf.net/python-collections-defaultdict/
    #subbar = []
    #subqual = []
    #subsam = []

    # Collect data for subset of barcodes
    ## "for" with "zip" is a good way to step through several lists in lockstep, instead of the traditional indexing
    ## zip takes two iterators, e.g.:
    ##      [0,1,2] and [a,b,c]
    ## and makes a new list like
    ##      [ (0,a), (1,b), (2,c) ]
    ## It's almost like Python's motto: when in doubt, use an iterator.
    for length, sam_row in zip(lengths, (barcodes, quals, samrows)):
    #for i in range(len(barcodes)):
    #    if lengths[i] == length:
    #        #if int(qual[i])>= 1:
            rows_by_length[length].append(sam_row)
            #subbar.append(barcodes[i])
            #subqual.append(qual[i])
            #subsam.append(samrows[i])

    # If subset of reads > 0 -> find duplicates
    # dict.items() is an iterator of key, value pairs
    samrows_dedup = []
    for length_subset, samrows_subset in rows_by_length.items():
        # Empty lists evaluate to False, so you can just use an if here. Whee! Python!
        if sam_rows:
        #if len(subbar) > 0:
            #samrows = identify_duplicates(subbar, subsam, subqual, samrows)
            samrows_dedup.append( identify_duplicates(samrows_subset) )
    #return samrows
    return samrows_dedup

# Identify duplicates and mark all but the one with the best alignment.
#def identify_duplicates(subbar, subsam, subqual, samrows):
def mark_duplicates(barcodes_list, lengths_list, quals_list, samrows_list):

    quals_dict = {}

    # sort by length, quality
    samrows_list.sort(samrows_list, key=lambda elt: (elt[8], elt[4])

    for barcode, qual, samrow in zip(barcodes_list, quals_list, samrows_list):
        ## "in" checks for the presence of barcode as a key in quals_dict. It is pythonic.
        ## Also numeric comparisons to None always favor of the number, so 4 > None is True
        if qual > quals_dict.get(barcode):
            quals_dict[barcode][ = qual

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
