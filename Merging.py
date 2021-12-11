#!/usr/bin/python


import os
import sys
from collections import Counter
import pandas as pd

input_path = sys.argv[sys.argv.index('-i')+1]
output_file = sys.argv[sys.argv.index('-o')+1]

liste = input_path.split(',')

premier = os.getcwd()

seconde = os.path.join(premier, output_file)

troisieme = os.path.join(premier, liste[0])

#quatrieme = os.listdir(troisieme)


#for i in quatrieme:
#	file = open(os.path.join(seconde, i), "w")

path_liste = []
final_liste = []


for l in liste:
	path_liste.append(os.path.join(premier, l))
	
for p in path_liste:
	file_liste = os.listdir(p)
	for f in file_liste:
		final_liste.append(p + "/" + f)



df = pd.DataFrame([f for f in final_liste], columns = ['Fullpath'])



file_names = df['Fullpath'].str.rsplit("/", 1, expand=True).rename(columns={0: 'Path', 1:'Filename'})


dfb = df.join(file_names)


for x in dfb['Filename']:
    paths = dfb[dfb['Filename'] == x]['Fullpath'] # Get list of fullpaths from unique filenames
    dfs = [path for path in paths] # Get list of dataframes from CSV file paths
    with open(os.path.join(seconde, x), 'w') as outfile:
    	for fname in dfs:
    		with open(fname) as infile:
    			for line in infile:
    				outfile.write(line)
counting = []
    
with open(os.path.join(seconde, "Manta_summary.txt"), "r") as outfile:
	for line in outfile:
		counting.append(line)

summary_output = open(os.path.join(seconde, "Concatenated_Manta_summary.txt"), "w")

Sums_of_counts = Counter(counting)
for key, value in Sums_of_counts.items():
    summary_output.write(key + '\t' + str(value) + '\n')






