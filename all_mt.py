 #!/usr/bin/env python

#ALL_MT - A pipeline developed to identify putative nuclear mitochondrial dna segments(NUMTS) in NGS contigs 
#This pipeline was developed by Carlos Henrique Aguiar Costa <carloscostaha@gmail.com>

import os
import sys
import argparse
import shutil
import subprocess
import pandas as pd
import numpy as np
import logging
import numt_modules as nm


#Needed Tools for Running ALL MT

programs_list = ['makeblastdb','blastn', 'seqtk', 'awk', 'samtools', 'bedtools', 'bowtie2', 'spades.py','megahit']

#Variables

BLASTN_1 = []
NUMTS_1 = []
NUMTS_2 = []
NUMTS_3 = []
NUMTS_4 = []
NUMTS_5 = []
numt_range = []
numt_list  = []
contigs_list = []
length = ''
unpaired = False

#Arguments

parser = argparse.ArgumentParser(description='ALLMT: A pipeline for identification of pseudogenes in mitocondrial data based in BLAST Identity, Assembly And Read Coverage',
	usage='%(prog)s [-h] -m [mitochondrial.fasta] -contigs [contigs.fasta] -assembler [spades OR megahit] -R1 [R1.fastq] -R2 [R2.fastq] -o [output_folder]', 
	epilog="Version: 1.0 Developed by: Carlos H. A. Costa - email: carloscostaha@gmail.com")

parser.add_argument('-m', '--mitochondrial', type = str, required = True, metavar = "", help = "Fasta file containing the first mitochondrial genome")
parser.add_argument('-contigs', '--contig_file', type = str, required = True,  metavar = "",help = "Fasta file containing the contigs from the metagenome assembly")
parser.add_argument('-assembler', type = str, required = True,  metavar = "",help="Please, inform which metagenome assembler generated your reads (e.g: 'spades', 'megahit')")
parser.add_argument('-R1', '--paired1', type = str, required = True,  metavar = "",help = "Fastq file containing the paired reads 5'")
parser.add_argument('-R2', '--paired2', type = str, required = True,  metavar = "",help = "Fastq file containing the paired reads 3'")
parser.add_argument('-U1', '--unpaired1', type = str, required = False,  metavar = "",help = "Fastq file containing the unpaired reads 5'")
parser.add_argument('-U2', '--unpaired2', type = str, required = False,  metavar = "",help = "Fastq file containing the unpaired reads 3'")
parser.add_argument('-db', '--blastdb', type = str, required = False,  metavar = "",help = "Path for DB containing the merged reads in fasta file for BLASTN")
parser.add_argument('-id','--identity', type = float, required = False,  metavar = "",help = "Set a float value (0-100) for identity threshold for BLASTN tool. Default: 95")
parser.add_argument('-o', '--output', type = str, required = True,  metavar = "",help = "Prefix For the OUTPUT Results Folder - output_RESULTS/")

args = parser.parse_args()

#Output Folder Creation and Logging Creation

folder = args.output #Output creation

if os.path.exists(folder) and os.path.isdir(folder):
	create = input("WARNING: Same folder name detected. Do you want to overwrite it?[Y/N]\n")
	if create in ['y','Y','Yes']:
		shutil.rmtree(folder)
		os.mkdir(folder)
	elif create in ['n','N','No']:
		folder = folder + "_new"
		os.mkdir(folder)
	else:
		print("Unavailable argument... exiting")
		exit()
else:
	os.mkdir(folder)

nm.logging_all_mt(folder)

# logging.basicConfig(
# 	format='%(asctime)s - %(message)s', 
# 	datefmt='%d-%b-%y %H:%M:%S',
# 	handlers = [logging.FileHandler('{0}/log.txt'.format(folder)),logging.StreamHandler()],
# 	level=logging.INFO)

logging.info('PROCESS HAS STARTED - ALL MT PIPELINE')

#Create a step for verification of installed softwares in the machine, if they're not installed, install.

for i in programs_list:
	nm.is_tool(i)

#Check identity value for BLASTN
if args.identity:
	identity_blastn = args.identity
else:
	identity_blastn = float(95)

#Runs BLASTN installed in the local machine against all the contigs from <-contigs>

cmd = "blastn -query {0} -subject {1} -outfmt \'6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand\' -max_hsps 1 -out {2}/blastn_mitochondrial.txt".format(args.mitochondrial, args.contig_file, folder)
os.system(cmd)

#Reading the file and insert headers
data = ("%s/blastn_mitochondrial.txt"  % (folder))
mitogenome_numts = pd.read_table('{0}'.format(data), header = None, sep = '\t', names = ["qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand"])

logging.info('STEP 1) Source Separation')
#print ("STEP 1) Source Separation")

mitochondrial_contigs = nm.mitogenome(mitogenome_numts, identity_blastn) 
numts_contigs_LVL1 = nm.numts_LVL1(mitogenome_numts, identity_blastn)

logging.info("STEP 2) Flanking region identification in 5'")
#print ("STEP 2) Flanking region identification in 5'") #If the sequece has both level it receives both LVL 2 and LVL 3 status, else, only one of them respectively - The LVL 2 stands for the beggining 5', and the LVL 3 for the end 3'

numts_contigs_LVL2 = nm.numts_LVL2(numts_contigs_LVL1)
#print (numts_contigs_LVL2)
logging.info("STEP 3) Flanking region identification in 3'")
#print ("STEP 3) Flanking region identification in 3'")

numts_contigs_LVL3 = nm.numts_LVL3(numts_contigs_LVL2)
logging.info("STEP 4) Contigs extraction and Remapping")
#print ("STEP 4) Contigs extraction and Remapping") #From the LVL 2 and LVL 3 contigs extracts the sequence ID, stores in temporary list, and maps every positive over the contigs pNUMTs
logging.info("\tSTEP 4.1) Creating Working Directory with Fasta Files and Bowtie2-DB")
#print ("\tSTEP 4.1) Creating Working Directory with Fasta Files and Bowtie2-DB")

#Creates a list which contains every positive pNUMTS so far
contigs_list = numts_contigs_LVL3["sseqid"].values.tolist()

#Creates a dir for each pNUMT based in the list
if args.assembler == "spades":
	logging.info('\t\t4.1.1)BlastDB with SPADES contigs')
	nm.dir_LVL_4_spades(contigs_list, args.contig_file,folder)
elif args.assembler == "megahit":
	logging.info('\t\t4.1.1)BlastDB with MEGAHIT contigs')
	nm.dir_LVL_4_megahit(contigs_list, args.contig_file,folder)
else:
	logging.error("No assembler was declared... exit")
	#print ("No assembler informed... exiting")
	exit()

#Creates a BLAST DB based in all the reads concatenated in fasta format, if not declared, the tool will use the declared reads;
if args.blastdb:
 	logging.info("\tSTEP 4.2) DB containing reads DECLARED:")
 	#print ("\tSTEP 4.2) DB containing reads DECLARED:")
 	blastn_file = args.blastdb
else:
 	blastn_file = nm.blastn_database(args.paired1, args.paired2, args.unpaired1, args.unpaired2,folder)

#From all the positive reads that matched with every single pNUMT, a list containing all of them is created 
logging.info("\tSTEP 4.3) Running BLAST against all the reads (in fasta) and creating a list with every unique read")
#print("\tSTEP 4.3) Running BLAST against all the reads (in fasta) and creating a list with every unique read")
nm.blastn_readlist(contigs_list, blastn_file,folder)	

#All the aligned read information is extracted from the original paired reads data.
logging.info("\tSTEP 4.4) Extracting the mapped reads - SEQTK")
#print ("\tSTEP 4.4) Extracting the mapped reads - SEQTK")
nm.blastn_reads(contigs_list, args.paired1, args.paired2, args.unpaired1, args.unpaired2,folder)

#Mapping all the previous aligned reads over their respective contigs (pNUMTS), return to the user a table with the coverage in each position
logging.info("\tSTEP 4.5) Mapping the Reads against the pNUMT - Bowtie2 and Samtools")
#print ("\tSTEP 4.5) Mapping the Reads against the pNUMT - Bowtie2 and Samtools")
if args.unpaired1 and args.unpaired2:
	unpaired = True

nm.bowtie_samtools_mapping(contigs_list, unpaired,folder)

logging.info("\tSTEP 4.6) Verification of GAPS inside the pNUMTs aligned regions")
#print ("\tSTEP 4.6) Verification of GAPS inside the pNUMTs aligned regions")

numts_LVL4 = nm.numts_LVL4(contigs_list,numts_contigs_LVL3,folder)
dataframe_list = pd.DataFrame(contigs_list) #This is a tric to put my contig list back to my dataframe, since in the previous step I've used to remove the data more easily and couldn't put it back... better investigate it.
numts_LVL4.insert(1,"sseqid",dataframe_list)

logging.info("STEP 5) Assembly of the the reads from every pNUMT sequence - ASSEMBLER")
#print ("STEP 5) Assembly of the the reads from every pNUMT sequence - ASSEMBLER")
logging.info("\tSTEP 5.1) Reassembling the pNUMTs")
#print ("\tSTEP 5.1) Reassembling the pNUMTs")

nm.numts_LVL5_reassembly(contigs_list, args.paired1, args.paired2, args.unpaired1, args.unpaired2, folder)

numts_LVL5 = nm.numts_LVL5(contigs_list,args.mitochondrial, numts_LVL4,folder)

logging.info("STEP 6) Final Level")
#print ("STEP 6) Final Level")

numts_max = nm.numts_max(numts_LVL5)

numts_max.to_csv("{0}/Final_Results.txt".format(folder), index = False, header = True, sep = '\t', float_format="%.2f")

print("\n----- Thank you for using ALL MT -----\n")

exit()

