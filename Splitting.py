#!/usr/bin/python2.7


import sys
import itertools

seqs_file = sys.argv[sys.argv.index('-i')+1]
output_file = sys.argv[sys.argv.index('-o')+1]

output = open(output_file + str(seqs_file[:-6]) + ".fasta", "w")

with open(seqs_file, 'r') as f:
    xlines = itertools.islice(f, 1, None, 4)
    for lines in xlines:
        output.write(lines)

        



		



