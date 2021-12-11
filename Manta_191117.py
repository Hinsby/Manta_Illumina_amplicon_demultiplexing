#!/usr/bin/python2.7

#For demultiplexing of dual-barcoded Illumina MiSeq data that has been merged prior with CASPER. Will accept multiple
#different primers in the same lane.

import re
import csv
from collections import Counter
import regex
import sys
import itertools
import os.path


#Usage is as follows: python Manta.py -i <path to fastq input file> -m <path to unix mapping file> -o <path to output folder>
#-c option for clipping barcodes and primer (T or F) -s integer for filtration sweeps, through 1 to 6

# python Manta_191117.py -i Lane4CasperMerged.fastqaa -m Lane4LinkSubsetMappingFile.txt -o Lane4LinkTest-cT/OutputAA -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqab -m Lane4LinkSubsetMappingFile.txt -o Lane4LinkTest-cT/OutputAB -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqac -m Lane4LinkSubsetMappingFile.txt -o Lane4LinkTest-cT/OutputAC -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqad -m Lane4LinkSubsetMappingFile.txt -o Lane4LinkTest-cT/OutputAD -c T -s 6


# python Manta_191117.py -i Lane4CasperMerged.fastqaa -m Lane4MappingFile9.24.18.txt -o 16SonlyAA -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqab -m Lane4MappingFile9.24.18.txt -o 16SonlyAB -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqac -m Lane4MappingFile9.24.18.txt -o 16SonlyAC -c T -s 6
# python Manta_191117.py -i Lane4CasperMerged.fastqad -m Lane4MappingFile9.24.18.txt -o 16SonlyAD -c T -s 6


seqs_file = sys.argv[sys.argv.index('-i')+1]
map_file = sys.argv[sys.argv.index('-m')+1]
output_file = sys.argv[sys.argv.index('-o')+1]
c = sys.argv[sys.argv.index('-c')+1]
s = sys.argv[sys.argv.index('-s')+1]

clipping = str(c)
sweeps = int(s)

	

print("                           **")
print("                          **")
print("                         * *")
print("                         *  *                    * *")
print("                         *    *               *")
print("                          *     *  * *      *")
print("                           *            * *")
print("                          *              *")
print("                         *                *")
print("                        *   *             *")
print("                              *           *")
print("                            *             *")
print("                            * * * *      *")
print("                                     *    *")
print("                                      *    *")
print("                                        *   *")
print("                                          *   *")
print("                                             * **")
print("Welcome to Manta, a demultiplexing tool that filter-feeds on merged")
print("DNA sequences, filtering by both barcode indices and primer pairs.")

#Begin by importing mapping file and creating lists of sample information.

def import_file(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator, skipinitialspace=True):
        if line:
            yield line

SampleIDs = [] #Column 1 of the mapping file: this corresponds to a unique sample name
Sample_forward_barcode = [] #Column 2 of the mapping file: this corresponds to a barcode of variable length (=< 10 bp)
Degen_forward_primer = [] #Column 3: this corresponds to a forward primer in the 5'-3' orientation. Degenerate bases in accord with standard IUPAC naming are acceptable.
Degen_reverse_primer = [] #Column 4: this corresponds to a reverse primer in the 5'-3' orientation.
Sample_reverse_barcode = [] #Column 5: this corresponds to a reverse barcode of variable length (=< 10 bp)

for row in import_file(map_file, '\t'):
    SampleIDs.append(row[0])
    Sample_forward_barcode.append(row[1])
    Degen_forward_primer.append(row[2])
    Degen_reverse_primer.append(row[3])
    Sample_reverse_barcode.append(row[4])

#Remove headers from the lists

SampleIDs.pop(0)
Sample_forward_barcode.pop(0)
Degen_forward_primer.pop(0)
Degen_reverse_primer.pop(0)
Sample_reverse_barcode.pop(0)

#Convert degenerate bases to regex format

iupac = {'[':'[', ']':']','A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
            'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
            'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}

Sample_forward_primer = []
Sample_reverse_primer = []

for curr_for_primer in Degen_forward_primer:
    Sample_forward_primer.append(''.join([iupac[symbol] for symbol in curr_for_primer]))
for curr_rev_primer in Degen_reverse_primer:
    Sample_reverse_primer.append(''.join([iupac[symbol] for symbol in curr_rev_primer]))


#Define reverse compliment function for list of barcodes

def reversecomplimentbarcode(barcode):
    for base in barcode:
        step_one = [b.replace("A", "X") for b in base]
        step_two = [t.replace("T", "Y") for t in step_one]
        step_three = [t.replace("G", "Z") for t in step_two]
        step_four = [t.replace("C", "W") for t in step_three]
        step_five = [t.replace("X", "T") for t in step_four]
        step_six = [t.replace("Y", "A") for t in step_five]
        step_seven = [t.replace("Z", "C") for t in step_six]
        step_eight = [t.replace("W", "G") for t in step_seven]
        Reverse = step_eight[::-1]
        Reverse = ''.join(Reverse)
        yield Reverse

Reverse_comp_forward_barcode = reversecomplimentbarcode(Sample_forward_barcode)
Reverse_comp_reverse_barcode = reversecomplimentbarcode(Sample_reverse_barcode)

#Create one object, SampleList, that has all SampleID, Sample_Barcode and Sample_primer info stored in it

SampleList = []

for i, j, k, l, m, n, o in zip(SampleIDs,Sample_forward_barcode,Sample_forward_primer, Sample_reverse_primer, Reverse_comp_reverse_barcode, Sample_reverse_barcode, Reverse_comp_forward_barcode):
    SampleList.append(">" + i + "\t" + j + "\t" + k + "\t" + l + "\t" + m + "\t" + n + "\t" + o)

#Import sequences as Seqs_liste.

Input_liste = []


with open(seqs_file, 'r') as f:
    for lines in f:
        Input_liste.append(lines)

Seqs_liste = Input_liste

Input_liste = []

#Define function to separate sequences: first based on presence of both forward and reverse barcodes within the first
#and last 10 bp of the read; second, based on presence of primers, using regex for primers due to degenerate bases.

#These lists must be constructed to store information during the various steps of filtering, whereby:
#Sorted_seqs is a liste for sequences that pass primer and barcode filtering.
#Summary_liste is an output of Sample names, which will be used to summarise total counts of identified sequences.
#Barcode_mismatches is an output of sequences which did not match the barcodes in the mapping file.

Sorted_seqs = []
First_sweep = []

seq_counter = 0


sample_counter = []

#Define a function to output the sequences into individual files that correspond to that SampleID name. 
#This allows for comparison of samples across sequencing lanes, although the author considers this to be 
#less than ideal if an experiment can be designed for comparative analyses within the same sequencing lane.

def Seq_outputting(IDs, Seqs):
	ID_counter = 0
	for i in IDs:
		file = open(os.path.join(output_file, i + ".fasta"), 'a+')
    		for s in Seqs:
        		seq_output = s.split('\t')[1] 
        		ID = s.split('\t')[0]
        		if i in ID:
					ID_counter += 1
					if clipping == "T":
						file.write(ID + "|" + str(ID_counter) + '\n' + seq_output[30:-30] + '\n')
						sample_counter.append(ID + ',')
					else:
    						file.write(ID + "|" + str(ID_counter) + '\n' + seq_output)
    						sample_counter.append(ID + ',')
				

summary_output = open(os.path.join(output_file, "Manta_summary.txt"), "w")

#Define a function for the output summary
def Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered):
	summary_output.write("Total sequences =" + '\t' + str(len(Seqs_liste)) + '\n')
	summary_output.write("Recovered sequences =" + '\t' + str(len(set(Sorted_seqs))) + '\n')
	Recovery_rate = str((float(len(set(Sorted_seqs))) / len(set(Seqs_liste)))*100)
	summary_output.write("Recovery rate =" + '\t' + str(Recovery_rate) + "(%)" + '\n')
	summary_output.write("Unassigned sequences =" + '\t' + str(len(Unrecovered)) + '\n')
	summary_output.write("Individual sample recovery:" + '\n')

#Define a function for binning the unassigned sequences
def Unassigned_output(Unrecovered):
	Barcode_mismatch_output = open(os.path.join(output_file, "Manta_unassigned_output.fasta"), "w")
	for mismatch in Unrecovered:
    		Barcode_mismatch_output.write(mismatch)
    		
#Define a function for counting the sequences per sample ID
def Counting_seqs(sample_counter):
	Sums_of_counts = Counter(sample_counter)		
	for key, value in Sums_of_counts.items():
    		summary_output.write(key + '\t' + str(value) + '\n')

Incorrect_bc_output = open(os.path.join(output_file, "Incorrect_barcodes.txt"), "w")
Incorrect_counter = []

#Define a function for identifying barcode issues in the unassigned bin, specifically if the barcodes don't match
#The mapping file due to human error
def Incorrect_barcodes(Unrecovered):
	for mismatch in Unrecovered:
		Inc_f_bc = mismatch[0:10]
		Inc_r_bc = mismatch[-10:]
		Incorrect_counter.append(Inc_f_bc + "_" + Inc_r_bc)
		

#And another function to count the incorrect barcodes		
def Counting_incorrect(Incorrect_counter):
	Stripped_line = [s.rstrip() for s in Incorrect_counter]
	Sums_of_counts = Counter(Stripped_line)
	for key, value in Sums_of_counts.items():
		Incorrect_bc_output.write(key + '\t' + str(value) + '\n')		


for sequence in Seqs_liste:
	if len(sequence) < 300: #Remove merged sequences less than 300 bp because they're likely not terribly wonderful
		Seqs_liste.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50] #This is where we hunt for our primer
        	front_bc_section = sequence[0:10] #And forward barcode
        	back_bc_section = sequence[-10:] #And reverse barcode
        	for_primer_match = regex.search('('+forpr+'){e<=3}', seq_section) #With regex, we allow up to 3 base mismatches
        	rev_primer_match = regex.search('('+revpr+'){e<=3}', seq_section)
        	forbc_match = regex.search('('+forbc+')', front_bc_section)
        	revbc_match = regex.search('('+revbc+')', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ ')', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + ')', back_bc_section)
        	if for_primer_match:
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		First_sweep.append(sequence) #Append sequences in 5'-3' orientation
        	elif rev_primer_match: #If in 3'-5' orientation, we reverse complement
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse) #Append other sequences in 5'-3'
                				First_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_first_pool = list(set(Seqs_liste)- set(First_sweep))


print("        **")
print("       *  *")
print("    ***     *")
print("  **  *     *")
print(" *     *  *")
print("        **")


print("First filtration stats")
print("Total input sequences:" + '\t' + str(len(Seqs_liste)))
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
First_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(First_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_first_pool)))
print("Recovered an extra:" + '\t' + str(len(First_sweep)))

First_sweep_output = open(os.path.join(output_file, "First_sweep.txt"), "w")
First_sweep_output.write("First filtration stats" + '\n')
First_sweep_output.write("Total input sequences:" + '\t' + str(len(Seqs_liste)) + '\n')
First_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
First_sweep_output.write("Misassigned sequences:" + '\t' + str(First_pool_error) + '\n')
First_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_first_pool)) + '\n')
First_sweep_output.write("Recovered an extra:" + '\t' + str(len(First_sweep)))



if sweeps == 1:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_first_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_first_pool)
	Incorrect_barcodes(Unrecovered_first_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass


Second_sweep = []

for sequence in Unrecovered_first_pool:
	if len(sequence) < 300:
		Unrecovered_first_pool.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50]
        	front_bc_section = sequence[1:9]
        	back_bc_section = sequence[-9:-1]
        	for_primer_match = regex.search('('+forpr+'){e<=3}', seq_section)
        	rev_primer_match = regex.search('('+revpr+'){e<=3}', seq_section)
        	forbc_match = regex.search('('+forbc+'){e<=1}', front_bc_section)
        	revbc_match = regex.search('('+revbc+'){e<=1}', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ '){e<=1}', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + '){e<=1}', back_bc_section)
        	if for_primer_match: #or rev_primer_match: #Filter via primer presence
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		Second_sweep.append(sequence)
        	elif rev_primer_match:
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse)
                				Second_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_second_pool = list(set(Unrecovered_first_pool)- set(Second_sweep))

print("    *      **")
print("     **   *  *")
print("       ***     *")
print("         *     *")
print("          *  *")
print("           **")


print("Second filtration stats")
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
Second_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(Second_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_second_pool)))
print("Recovered an extra:" + '\t' + str(len(Second_sweep)))

Second_sweep_output = open(os.path.join(output_file, "Second_sweep.txt"), "w")
Second_sweep_output.write("Second filtration stats" + '\n')
Second_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
Second_sweep_output.write("Misassigned sequences:" + '\t' + str(Second_pool_error) + '\n')
Second_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_second_pool)) + '\n')
Second_sweep_output.write("Recovered an extra:" + '\t' + str(len(Second_sweep)))

if sweeps == 2:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_second_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_second_pool)
	Incorrect_barcodes(Unrecovered_second_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass


Third_sweep = []

Unrecovered_first_pool = []

for sequence in Unrecovered_second_pool:
	if len(sequence) < 300:
		Unrecovered_second_pool.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50]
        	front_bc_section = sequence[1:10]
        	back_bc_section = sequence[-10:-1]
        	for_primer_match = regex.search('('+forpr+'){e<=3}', seq_section)
        	rev_primer_match = regex.search('('+revpr+'){e<=3}', seq_section)
        	forbc_match = regex.search('('+forbc+'){e<=1}', front_bc_section)
        	revbc_match = regex.search('('+revbc+'){e<=1}', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ '){e<=1}', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + '){e<=1}', back_bc_section)
        	if for_primer_match: #or rev_primer_match: #Filter via primer presence
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		Third_sweep.append(sequence)
        	elif rev_primer_match:
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse)
                				Third_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_third_pool = list(set(Unrecovered_second_pool)- set(Third_sweep))

print("              **")
print("             *  *")
print("          ***     *")
print("        **  *     *")
print("       *     *  *")
print("              **")


print("Third filtration stats")
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
Third_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(Third_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_third_pool)))
print("Recovered an extra:" + '\t' + str(len(Third_sweep)))

Third_sweep_output = open(os.path.join(output_file, "Third_sweep.txt"), "w")
Third_sweep_output.write("Third filtration stats" + '\n')
Third_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
Third_sweep_output.write("Misassigned sequences:" + '\t' + str(Third_pool_error) + '\n')
Third_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_third_pool)) + '\n')
Third_sweep_output.write("Recovered an extra:" + '\t' + str(len(Third_sweep)))

if sweeps == 3:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_third_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_third_pool)
	Incorrect_barcodes(Unrecovered_third_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass

Fourth_sweep = []

Unrecovered_second_pool = []

for sequence in Unrecovered_third_pool:
	if len(sequence) < 300:
		Unrecovered_third_pool.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50]
        	front_bc_section = sequence[1:12]
        	back_bc_section = sequence[-12:-1]
        	for_primer_match = regex.search('('+forpr+'){e<=3}', seq_section)
        	rev_primer_match = regex.search('('+revpr+'){e<=3}', seq_section)
        	forbc_match = regex.search('('+forbc+'){e<=1}', front_bc_section)
        	revbc_match = regex.search('('+revbc+'){e<=1}', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ '){e<=1}', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + '){e<=1}', back_bc_section)
        	if for_primer_match: #or rev_primer_match: #Filter via primer presence
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		Fourth_sweep.append(sequence)
        	elif rev_primer_match:
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse)
                				Fourth_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_fourth_pool = list(set(Unrecovered_third_pool)- set(Fourth_sweep))

print("          *      **")
print("           **   *  *")
print("             ***     *")
print("               *     *")
print("                *  *")
print("                 **")


print("Fourth filtration stats")
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
Fourth_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(Fourth_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_fourth_pool)))
print("Recovered an extra:" + '\t' + str(len(Fourth_sweep)))

Fourth_sweep_output = open(os.path.join(output_file, "Fourth_sweep.txt"), "w")
Fourth_sweep_output.write("Fourth filtration stats" + '\n')
Fourth_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
Fourth_sweep_output.write("Misassigned sequences:" + '\t' + str(Fourth_pool_error) + '\n')
Fourth_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_fourth_pool)) + '\n')
Fourth_sweep_output.write("Recovered an extra:" + '\t' + str(len(Fourth_sweep)))

if sweeps == 4:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_fourth_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_fourth_pool)
	Incorrect_barcodes(Unrecovered_fourth_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass

Fifth_sweep = []

Unrecovered_third_pool = []

for sequence in Unrecovered_fourth_pool:
	if len(sequence) < 300:
		Unrecovered_fourth_pool.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50]
        	front_bc_section = sequence[0:12]
        	back_bc_section = sequence[-12:]
        	for_primer_match = regex.search('('+forpr+'){e<=3}', seq_section)
        	rev_primer_match = regex.search('('+revpr+'){e<=3}', seq_section)
        	forbc_match = regex.search('('+forbc+'){e<=1}', front_bc_section)
        	revbc_match = regex.search('('+revbc+'){e<=1}', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ '){e<=1}', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + '){e<=1}', back_bc_section)
        	if for_primer_match: #or rev_primer_match: #Filter via primer presence
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		Fifth_sweep.append(sequence)
        	elif rev_primer_match:
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse)
                				Fifth_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_fifth_pool = list(set(Unrecovered_fourth_pool)- set(Fifth_sweep))

print("                    **")
print("                   *  *")
print("                ***     *")
print("              **  *     *")
print("             *     *  *")
print("                    **")


print("Fifth filtration stats")
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
Fifth_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(Fifth_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_fifth_pool)))
print("Recovered an extra:" + '\t' + str(len(Fifth_sweep)))

Fifth_sweep_output = open(os.path.join(output_file, "Fifth_sweep.txt"), "w")
Fifth_sweep_output.write("Fifth filtration stats" + '\n')
Fifth_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
Fifth_sweep_output.write("Misassigned sequences:" + '\t' + str(Fifth_pool_error) + '\n')
Fifth_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_fifth_pool)) + '\n')
Fifth_sweep_output.write("Recovered an extra:" + '\t' + str(len(Fifth_sweep)))

if sweeps == 5:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_fifth_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_fifth_pool)
	Incorrect_barcodes(Unrecovered_fifth_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass

Sixth_sweep = []

Unrecovered_fourth_pool = []

for sequence in Unrecovered_fifth_pool:
	if len(sequence) < 300:
		Unrecovered_fifth_pool.remove(sequence)
    	for sample in SampleList:
        	ID, forpr, revpr, forbc, revbc, R_revbc, R_forbc = sample.split('\t')[0],sample.split('\t')[2], sample.split('\t')[3], sample.split('\t')[1], sample.split('\t')[4], sample.split('\t')[5], sample.split('\t')[6]
        	seq_section = sequence[0:50]
        	front_bc_section = sequence[0:12]
        	back_bc_section = sequence[-12:]
        	for_primer_match = regex.search('('+forpr+'){e<=5}', seq_section)
        	rev_primer_match = regex.search('('+revpr+'){e<=5}', seq_section)
        	forbc_match = regex.search('('+forbc+'){e<=1}', front_bc_section)
        	revbc_match = regex.search('('+revbc+'){e<=1}', back_bc_section)
        	R_forbc_match = regex.search('(' +R_revbc+ '){e<=1}', front_bc_section)
        	R_revbc_match = regex.search('(' + R_forbc + '){e<=1}', back_bc_section)
        	if for_primer_match: #or rev_primer_match: #Filter via primer presence
            		if forbc_match and revbc_match: #Begin filtering by combinations of barcode location
                		Sorted_seqs.append(ID + '\t' + sequence)
                		Sixth_sweep.append(sequence)
        	elif rev_primer_match:
            			if R_forbc_match and R_revbc_match:
        					step_one = [base.replace("A", "X") for base in sequence]
        					step_two = [t.replace("T", "Y") for t in step_one]
        					step_three = [t.replace("G", "Z") for t in step_two]
        					step_four = [t.replace("C", "W") for t in step_three]
        					step_five = [t.replace("X", "T") for t in step_four]
        					step_six = [t.replace("Y", "A") for t in step_five]
        					step_seven = [t.replace("Z", "C") for t in step_six]
        					step_eight = [t.replace("W", "G") for t in step_seven]
        					Reverse = step_eight[::-1]
        					Reverse = ''.join(Reverse)
                				Sorted_seqs.append(ID + '\t' + Reverse)
                				Sixth_sweep.append(sequence)
        	else:
        		pass
    	else:
        	pass

Unrecovered_sixth_pool = list(set(Unrecovered_fifth_pool)- set(Sixth_sweep))

print("                *      **")
print("                 **   *  *")
print("                   ***     *")
print("                     *     *")
print("                      *  *")
print("                       **")


print("Sixth filtration stats")
print("Recovered sequences:" + '\t' + str(len(Sorted_seqs)))
Sixth_pool_error = len(Sorted_seqs) - len(set(Sorted_seqs))
print("Misassigned sequences:" + '\t' + str(Sixth_pool_error))
print("Unassigned sequences:" + '\t' + str(len(Unrecovered_sixth_pool)))
print("Recovered an extra:" + '\t' + str(len(Sixth_sweep)))

print("                          **")
print("                         *  *")
print("                      ***     *")
print("                    **  *     *")
print("                   *     *  *")
print("                          **")

Sixth_sweep_output = open(os.path.join(output_file, "Sixth_sweep.txt"), "w")
Sixth_sweep_output.write("Sixth filtration stats" + '\n')
Sixth_sweep_output.write("Recovered sequences:" + '\t' + str(len(Sorted_seqs)) + '\n')
Sixth_sweep_output.write("Misassigned sequences:" + '\t' + str(Sixth_pool_error) + '\n')
Sixth_sweep_output.write("Unassigned sequences:" + '\t' + str(len(Unrecovered_sixth_pool)) + '\n')
Sixth_sweep_output.write("Recovered an extra:" + '\t' + str(len(Sixth_sweep)))

if sweeps == 6:
	print("Collating outputs.")
	Seq_outputting(SampleIDs, Sorted_seqs)
	Summary_outputting(Seqs_liste, Sorted_seqs, Unrecovered_sixth_pool)
	Counting_seqs(sample_counter)
	Unassigned_output(Unrecovered_sixth_pool)
	Incorrect_barcodes(Unrecovered_sixth_pool)
	Counting_incorrect(Incorrect_counter)
	sys.exit()
else:
	pass

