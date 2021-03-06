#!/bin/bash
#Load in python 3.6.4, so you can run virtualenv
module load python/3.6.4
#load in the virtual enviorment you need to run Manta
source ~/virtual-envs/MANTA/bin/activate
#load in casper for splitting
module load casper/0.8.2
#Set global variable mapfile equal to the 3rd argument passed through shell script
export mapfile=$3

#IMPORTANT: MAPPING FILE MUST HAVE TAB-DELIMITED HEADER/COLUMN STRUCTURE AS SHOWN EXACTLY BELOW:
#SampleID	FwdBarcode	FwdTargetPrimer	RevTargetPrimer	RevBarcodeSeq

#use casper to generate a paired file, and remove reads that are unpaired
casper $1 $2 -t 28 -l
#convert casper fastq file too a fasta file
python Splitting.py -i casper.fastq -o ./
#split the fasta file into 28 equally sized files to speed up Manta processing time
split -l $(( $( wc -l < casper.fasta ) / 28 + 1 )) casper.fasta Manta_Input #MR changed from 14 to 28 so that Manta splitted inputs match the 28 output directories made below
#make directorys for all output files generated by Manta, as well as a final directory to stored the merged files
mkdir OutputAA
mkdir OutputAB
mkdir OutputAC
mkdir OutputAD
mkdir OutputAE
mkdir OutputAF
mkdir OutputAG
mkdir OutputAH
mkdir OutputAI
mkdir OutputAJ
mkdir OutputAK
mkdir OutputAL
mkdir OutputAM
mkdir OutputAN
mkdir OutputAO
mkdir OutputAP
mkdir OutputAQ
mkdir OutputAR
mkdir OutputAS
mkdir OutputAT
mkdir OutputAU
mkdir OutputAV
mkdir OutputAW
mkdir OutputAX
mkdir OutputAY
mkdir OutputAZ
mkdir OutputBA
mkdir OutputBB
mkdir MergedReads
#using the command parallel, and the file ALLjobs2run, manta will utilize all 28 threads requested and process each of the 14 splitted files at the same time
parallel --jobs 28 < ALLjobs2run
#merge the output files from all of the 14 individual Manta runs 
python Merging.py -i OutputAA,OutputAB,OutputAC,OutputAD,OutputAE,OutputAF,OutputAG,OutputAH,OutputAI,OutputAJ,OutputAK,OutputAL,OutputAM,OutputAN,OutputAO,OutputAP,OutputAQ,OutputAR,OutputAS,OutputAT,OutputAU,OutputAV,OutputAW,OutputAX,OutputAY,OutputAZ,OutputBA,OutputBB -o MergedReads
