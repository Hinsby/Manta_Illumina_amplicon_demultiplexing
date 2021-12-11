# Manta_Illumina_amplicon_demultiplexing
all the scripts needed to run the in-house demultiplexing scripts (collectively termed MANTA)

A how to guide for MANTA usage
Primary Author: Damien Finn, Ph.D. and Hinsby Cadillo-Quiroz
PI: Hinsby Cadillo-Quiroz
Funded by: Cadillo Lab at Arizona State University
Contributing Author: , Patrick Browne (build first version updaed by D Finn), Michael Pavia
Code collated to be presented to the public by: Mark Reynolds


Step 1: [assuming to run on computer cluster instead of personal computer] Create a virtual environment and install specific packages within this environment.
#Note MANTA was designed to be multi-threaded on 28 computer cores to increase performance speed. 

$mkdir ~/virtual-envs/
$virtualenv -p python2.7 ~/virtual-envs/MANTA
$source ~/virtual-envs/MANTA/bin/activate
$module load gcc/9.2.0 #used specifically for ASU's computer cluster; will not be applicable to new user outside of organization.
$pip install pandas
$pip install regex

#Running on a personal computer is discouraged since the script is designed to use 28 computational cores and requires this infrastructure.

Step 2: Transfer the following scripts used for processing and put them in the same folder that contains the data to be demultiplexed.

Alljobs2run
Manta_191117.py
Merging.py
NGS_Read_Processing_entire_for_28paralleljobs.sh
Splitting.py
Mapping_File.txt (this contains the barcode and sample information and changes for the paired end, multiplexed sequencing lane you which to demultiplex)

Step 3: Ensure scripts are executable.

$chmod 777 ALLjobs2run
$chmod 777 Manta_191117.py
$chmod 777 Merging.py
$chmod 777 NGS_Read_Processing_entire_for_28paralleljobs.sh
$chmod 777 Splitting.py

Step 4: Run MANTA.

$cd /directory/that/contains/everything
$bash NGS_Read_Processing_entire_for_28paralleljobs.sh Forward_Reads.fastq Reverse_Reads.fastq Mapping_file.txt

#Note it is wise to plan for a Illumina V3 2x300 lane with 150-200 samples for the run to take 48-72 hours to complete.

Appendix:

Explanation of NGS_Read_Processing_entire_for_28paralleljobs.sh:

The file packages the entire demultiplexing pipeline into and easy to run shell that runs each step in order.

1.	Loads in python version 3.6.4, where the program virtualenv is installed. 
2.	It will then load in the virtual environment MANTA (that you build previously) to ensure the correct versions of MANTAs dependencies are accessible to the computer
3.	Loads in, post-CASPER, per each of the 28 inputs and run using 28 threads (this makes it go faster), merges the pair ends from each input. This program matches and pairs reads from the forward and reverse fastq files. 
#See more here and the associated citation: http://best.snu.ac.kr/casper/index.php?name=manual
4.	It will then use the splitting.py script to convert the CASPER fastq output into a raw sequence file where each line is associated to an individual read. 
5.	Using the command split, this file is then separated into 28 different files. This step is crucial to increasing the speed of running MANTA. 
6.	Then 28 Output folders are made for MANTA outputs.
7.	The command parallel allows you to take advantage of the 28 threads that you reserved in you bash script. Written in the ALLjobs2run file contains 28 commands for running MANTA on all 28 subsets of the original data. This means that that the computer will simultaneously run MANTA on all these subsets and drastically decrease the time it takes to process the data. 
8.	Finally the Merging.py command will merge the outputs from each of the 28 MANTA runs and output the combined data into one folder names MergedReads. 


This piece is release as it is, for any bug or not expected behaviour please reach authors
