import os
import sys
import argparse
import shutil
import subprocess
import logging
from typing_extensions import Concatenate
import numpy as np
import pandas as pd

#Logging

def logging_all_mt(directory):
	logging.basicConfig(
		format='%(asctime)s - %(message)s', 
		datefmt='%d-%b-%y %H:%M:%S',
		handlers = [logging.FileHandler('{0}/log.txt'.format(directory)),logging.StreamHandler()],
		level=logging.INFO)

#Modules and Functions

#Function to read in the BLAST File the mitogenomes, separates what is a pnumt from a mitochondrial genome.

#Checking if tools are installed in the system and inside PATH variable

def is_tool(name):
	from whichcraft import which
	if which(name) is None:
		logging.error("{} is not installed, or found in the PATH variable".format(name))
		#print ("{} is not installed, or found in the PATH variable".format(name))
		exit()

#Step 1) Separation between putative mitochondrial contigs and pNUMTS, based upon identity #Pandas
def mitogenome(file1, identity): #If values from BLAST record are higher than set identity (or 0.95) will call it as mitochondrial contig
	mitochondrial_contigs  = file1.loc[file1['pident'] >= identity]
	mitochondrial_contigs = mitochondrial_contigs.copy()
	mitochondrial_contigs["LVL1"] = 0
	return file1

def numts_LVL1(file1, identity): #If values from BLAST record are lower than set identity (or 0.95) will call it as pNUMT
	numts_copied = file1.loc[file1['pident'] < identity]
	file1 = numts_copied.copy()
	file1["LVL1"] = 1
	return file1

#Step 2) Identification of nuclear flanking regions in 5' in each contig

def numts_LVL2(numts_LVL1, flank):
	try:
		pos = []
		neg = []
		numts_pos = numts_LVL1.loc[numts_LVL1['sstrand'].str.match('plus')]
		pos = numts_pos.copy()
		pos["LVL2"] = pos["sstart"].apply(lambda x: '1' if x >= flank else '0')
		neg = numts_LVL1.loc[numts_LVL1['sstrand'].str.match('minus')]
		neg = neg.copy()
		neg["LVL2"]= neg["send"].apply(lambda x: '1' if x >= flank else '0')
		concatenate = [pos, neg]
		results = pd.concat(concatenate)
		return results
	except ValueError:
		logging.error("\tError: There is no sequences under the identity (-id) threshold, please set a higher value")
		#print ("\tError: There is no sequences under the identity (-id) threshold, please set a higher value")
		exit()

#Step 3) Identification of nuclear flanking regions in 3' in each contig

#The number are not being properly reported in LVL3
#For plus strand --- if the length of contig minus 100 is greater than final position of the alignment, print "1", else "0"
#For minus strand --- if the length of contig minus 100 is greater than initial position of the alignment, print "1", else "0"

def numts_LVL3(numts_LVL2, flank):
	try:
		pos = []
		neg = []
		pos = numts_LVL2.loc[numts_LVL2["sstrand"].str.match('plus')].copy()
		pos.loc[pos["slen"] - pos["send"] >= flank, 'LVL3'] = "1"
		pos.loc[pos["slen"] - pos["send"] < flank, 'LVL3'] = "0"
		neg = numts_LVL2.loc[numts_LVL2["sstrand"].str.match('minus')].copy()
		neg.loc[neg["slen"] - neg["sstart"] >= flank, 'LVL3'] = "1"
		neg.loc[neg["slen"] - neg["sstart"] < flank, 'LVL3'] = "0"
		concat = [pos, neg]
		results = pd.concat(concat)
		return results	
	except ValueError:
		logging.error("\tError: There is no sequences under the identity (-id) threshold, please set a higher value")
		#print ("\tError: There is no sequences under the identity (-id) threshold, please set a higher value")
		exit()
#Step 4.1) Dir creation of every single pNUMT based in the previous results.

def dir_LVL_4_spades(numts_list, contig_file,folder): #creating subdirectory for each positive contig from SPADES assembly
	for i in numts_list:
		contig_dir = "{0}/{1}".format(folder,i)
		if os.path.exists(contig_dir) and os.path.isdir(contig_dir): #Looks if the directory file exist for every contig in the main dir.
			shutil.rmtree(contig_dir)	#Deletes every dir and it's previous content
			os.mkdir(contig_dir)	#Recreates the dir
			cmd = "awk -vRS='>' '{if(/%s/){print \">\" $0}}' %s > %s/%s/%s.fasta" % (i,contig_file,folder,i,i)
			os.system(cmd)
			os.system("mkdir {0}/{1}/bowtie2".format(folder, i))
			bowtie2 = "bowtie2-build %s/%s/%s.fasta %s/%s/bowtie2/%s_db -q 2> /dev/null" % (folder,i,i,folder,i,i)
			os.system(bowtie2)
		else:
			os.mkdir(contig_dir)
			cmd = "awk -vRS='>' '{if(/%s/){print \">\" $0}}' %s > %s/%s/%s.fasta" % (i,contig_file,folder,i,i)
			os.system(cmd)
			os.system("mkdir {0}/{1}/bowtie2".format(folder, i))			
			bowtie2 = "bowtie2-build %s/%s/%s.fasta %s/%s/bowtie2/%s_db -q 2> /dev/null"% (folder,i,i,folder,i,i)
			os.system(bowtie2)

def dir_LVL_4_megahit(numts_list, contig_file,folder): #creating subdirectory for each positive contig from MEGAHIT assembly
	for i in numts_list:
		contig_dir = "{0}/{1}".format(folder,i)
		if os.path.exists(contig_dir) and os.path.isdir(contig_dir): #Looks if the directory file exist for every contig in the main dir.
			shutil.rmtree(contig_dir)	#Deletes every dir and it's previous content
			os.mkdir(contig_dir)	#Recreates the dir
			cmd = "awk -vRS='>' '{if(/%s /){print \">\" $0}}' %s > %s/%s/%s.fasta" % (i,contig_file,folder,i,i)
			os.system(cmd)
			os.system("mkdir {0}/{1}/bowtie2".format(folder, i))
			bowtie2 = "bowtie2-build %s/%s/%s.fasta %s/%s/bowtie2/%s_db -q 2> /dev/null" % (folder,i,i,folder,i,i)
			os.system(bowtie2)
		else:
			os.mkdir(contig_dir)
			cmd = "awk -vRS='>' '{if(/%s /){print \">\" $0}}' %s > %s/%s/%s.fasta" % (i,contig_file,folder,i,i)
			os.system(cmd)
			os.system("mkdir {0}/{1}/bowtie2".format(folder, i))
			bowtie2 = "bowtie2-build %s/%s/%s.fasta %s/%s/bowtie2/%s_db -q 2> /dev/null"% (folder,i,i,folder,i,i)
			os.system(bowtie2)			

#Try do reduce with BioSeq::IO - Building the database with all the reads
def blastn_database(paired1, paired2, unpaired1, unpaired2,folder):
	logging.info("\tSTEP 4.2) DB containing reads NOT declared: Creating database as DB for further analysis")
	os.system("mkdir {0}/database_folder".format(folder))
	if unpaired1 and unpaired2:
		logging.info("\t\tSTEP 4.2.1)Creating unique temporary FASTQ File")
		total_reads = "cat {0} {1} {2} {3} > {4}/database_folder/total_reads.fastq".format(paired1, paired2, unpaired1, unpaired2, folder)
		os.system(total_reads)
	else:
		logging.info("\t\tSTEP 4.2.1)Creating unique temporary FASTQ File")	
		total_reads = "cat {0} {1} > {2}/database_folder/total_reads.fastq".format(paired1, paired2,folder)
		os.system(total_reads)
	logging.info("\t\tSTEP 4.2.2)Converting FASTQ File into FASTA File")
	seqtk = "seqtk seq -A {0}/database_folder/total_reads.fastq > {0}/database_folder/total_reads.fasta".format(folder) 
	os.system(seqtk)
	del_fastq = "rm {0}/database_folder/total_reads.fastq".format(folder) #deletes fastq file
	os.system(del_fastq)	
	logging.info("\t\tSTEP 4.2.3)Creating reads BLASTDB")
	blastn_db = "makeblastdb -in {0}/database_folder/total_reads.fasta -dbtype nucl -input_type fasta -out {0}/database_folder/DB 2> /dev/null".format(folder)
	os.system(blastn_db)		
	blastn_file = "{0}/database_folder/DB".format(folder)
	return blastn_file

#A blast of the putative numt - query, and all the reads -subject. From this results, a list containing only the single ID of the reads will be created.
def blastn_readlist(numts_list, blast_db,folder):
	for x in numts_list:
		blastn = "blastn -query %s/%s/%s.fasta -db %s -num_alignments 80000000 -outfmt '6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' -perc_identity 100 | awk '{if($4 == $6 && $6 >=  50){print $0}}' > %s/%s/%s_blastn_results.txt" % (folder,x,x,blast_db,folder,x,x)
		mapped_reads = "cut -f 3 %s/%s/%s_blastn_results.txt | sort | uniq > %s/%s/%s_read_list.txt" % (folder,x,x,folder,x,x)
		os.system(blastn)
		os.system(mapped_reads)

#This function takes all the paired, or not, reads, and extract only the ones whose hits where higher than 100% in the previous step, and produces, for each one of them, a file containing only those reads for further steps.
def blastn_reads(numts_list,paired1, paired2, unpaired1, unpaired2,folder):
	for i in numts_list:
	 	if unpaired1 and unpaired2:
	 		reads_paired_1 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R1P.fastq".format(paired1,folder,i)
	 		reads_paired_2 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R2P.fastq".format(paired2,folder,i)
	 		reads_unpaired_1 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R1U.fastq".format(unpaired1,folder,i)
	 		reads_unpaired_2 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R2U.fastq".format(unpaired2,folder,i)
	 		os.system(reads_paired_1)
	 		os.system(reads_paired_2)
	 		os.system(reads_unpaired_1)
	 		os.system(reads_unpaired_2)

	 	else:
	 		reads_paired_1 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R1P.fastq".format(paired1,folder,i)
	 		reads_paired_2 = "seqtk subseq {0} {1}/{2}/{2}_read_list.txt > {1}/{2}/{2}_R2P.fastq".format(paired2,folder,i)
	 		os.system(reads_paired_1)
	 		os.system(reads_paired_2)

#The reads are mapped over their respective pNUMTs and raw coverage is produced in order to look for putative gaps
def bowtie_samtools_mapping(numts_list, unpaired,folder):
	for seqs in numts_list:
		if unpaired == True:
			bowtie2_mapping = "bowtie2 -x {0}/{1}/bowtie2/{1}_db -1 {0}/{1}/{1}_R1P.fastq -2 {0}/{1}/{1}_R2P.fastq -U {0}/{1}/{1}_R1U.fastq -U {0}/{1}/{1}_R2U.fastq  --quiet --no-unal -a | samtools view -bS@ 32 - > {0}/{1}/{1}.bam 2> /dev/null".format(folder,seqs)
		else:
			bowtie2_mapping = "bowtie2 -x {0}/{1}/bowtie2/{1}_db -1 {0}/{1}/{1}_R1P.fastq -2 {0}/{1}/{1}_R2P.fastq --quiet --no-unal -a | samtools view -bS@ 32 - > {0}/{1}/{1}.bam 2> /dev/null".format(folder, seqs)
		samtools = "samtools sort -n {0}/{1}/{1}.bam > {0}/{1}/{1}_sort.bam 2> /dev/null".format(folder,seqs)
		bedtools = "bedtools genomecov -ibam {0}/{1}/{1}_sort.bam -d > {0}/{1}/{1}_pNUMT_coverage.txt".format(folder, seqs)
		os.system(bowtie2_mapping)
		os.system(samtools)
		os.system(bedtools)

def numts_LVL4(numts_list, numts_LVL3,folder):
	indexed = numts_LVL3.set_index("sseqid")
	results = []
	for i in numts_list:
		coverage = pd.read_table('{0}/{1}/{1}_pNUMT_coverage.txt'.format(folder,i), names = ["sequence", "position", "depth_coverage"])
		gapped_region = coverage[coverage["depth_coverage"] == 0] #Creates a copy of the list, with every gap occurence in the file
		gaps = gapped_region["position"].values.tolist()
		is_gapped = 0
		if gaps:
			gap = indexed.loc[i].copy()
			if  gap['sstrand'] == 'plus':
				start = gap["sstart"]
				end = gap["send"]
				for i in gaps:
					if start <= i <= end:
						is_gapped += 1
				if is_gapped >= 1:
					gap['LVL4'] = 0
				else:
					gap['LVL4'] = 0.5
			else:
				start = gap["send"]
				end = gap["sstart"]
				for i in gaps:
					if start <= i <= end:
						is_gapped += 1
				if is_gapped >= 1:
					gap['LVL4'] = 0
				else:
					gap['LVL4'] = 0.5
			results.append(gap)
		else:
			full = indexed.loc[i].copy()
			full["LVL4"] = 1
			results.append(full)
	NUMT_LVL4 = pd.concat(results, axis = 1, ignore_index = True)
	NUMT_LVL4 = NUMT_LVL4.T
	return NUMT_LVL4

#Starts to reassembly the reads, at first sight will use only MEGAHIT
#Some feature may change overtime: ASSEMBLY TYPE - I want to implement more than just one method, and also k-mer size control for each approach.
#debugging
def numts_LVL5_reassembly(numts_list, paired1, paired2, unpaired1,unpaired2,folder):
	if unpaired1 and unpaired2: #Process the paired and unpaired reads - R1;R2;U1;U2
		for sequence in numts_list:
			contig_assembly = "megahit -1 {0}/{1}/{1}_R1P.fastq -2 {0}/{1}/{1}_R2P.fastq -r {0}/{1}/{1}_R1U.fastq -r {0}/{1}/{1}_R2U.fastq --k-list 21,33,55,77,99,121 -o {0}/{1}/megahit --out-prefix {1} 2> /dev/null".format(folder,sequence)
			os.system(contig_assembly)
	else: #Process only the PAIRED reads - R1;R2
		for sequence in numts_list:
			contig_assembly = "megahit -1 {0}/{1}/{1}_R1P.fastq -2 {0}/{1}/{1}_R2P.fastq --k-list 21,33,55,77,99,121 -o {0}/{1}/megahit --out-prefix {1} 2> /dev/null".format(folder,sequence)
			os.system(contig_assembly)


def numts_LVL5(numts_list, mitochondrial_genome, numts_LVL4, folder):
	LVL5 =[]
	for sequence in numts_list: #Performs a blast of the pNUMT contigs against the newly produced contigs from the reads that belongs to it. Filtering only the best hist for each possible contig. The ideal scenario is to have a unique contig.
		blast_pnumts = "blastn -query {0}/{1}/{1}.fasta -subject {0}/{1}/megahit/{1}.contigs.fa -outfmt '6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand score' -max_hsps 1 > {0}/{1}/{1}_blast_score.txt 2> /dev/null".format(folder, sequence)
		os.system(blast_pnumts)
		#Parse the pNUMT against the assembled contig.
		blast_table = pd.read_table("{0}/{1}/{1}_blast_score.txt".format(folder,sequence), names = ["qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand", "score"])
		if blast_table.empty: #If there's no sequence inside the file
			LVL5.append(0)
		else:
			if len(blast_table.index) == 1:
				sequence_set = numts_LVL4[numts_LVL4["sseqid"] == sequence]
				qstart = sequence_set["qstart"].to_string(index = False) # query start
				qend = sequence_set["qend"].to_string(index = False) # query end
				contig_blast = 'blastn -query {0} -subject {1}/{2}/megahit/{2}.contigs.fa -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" > {1}/{2}/blast_mitogenome_contigs.txt 2> /dev/null'.format(mitochondrial_genome, folder, sequence)				
				os.system(contig_blast)
				blast_re = pd.read_table("{0}/{1}/blast_mitogenome_contigs.txt".format(folder,sequence), header = None, sep = '\t', names = ["qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand"])
				cstart = blast_re["qstart"].to_string(index = False) #contig start
				cend = blast_re["qend"].to_string(index = False) # contig end
				if cstart <= qstart:
					if cend >= qend:
						LVL5.append(1)
					else:
						LVL5.append(0)
				else:
					LVL5.append(0)	
			elif len(blast_table.index) > 1: #In the case of multiple contigs, this can occur if: 1) The contig is broken - In that case in STEP 4 it'll receive a value of 0.5; 2) If there's more than 1 mitochondrial genome in the sample, which might not align with the given reference.
				contigs_list = blast_table["sseqid"].values.tolist()
				contig_blast = 'blastn -query {0} -subject {1}/{2}/megahit/{2}.contigs.fa -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" > {1}/{2}/blast_mitogenome_contigs.txt 2> /dev/null'.format(mitochondrial_genome, folder, sequence)				
				os.system(contig_blast) #A blast of the multiple contigs against the mitochondrial genome, in order to identify if the pNUMT region was retrieved (qstart and qend). If positive, the process ends here, and the sequence receives a value of 0.5, if not found, the sequence receives de value of 0, assuming it could not be assembled.
				blast_re = pd.read_table("{0}/{1}/blast_mitogenome_contigs.txt".format(folder, sequence), header = None, sep = '\t', names = ["qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand"])
				sequence_set = numts_LVL4[numts_LVL4["sseqid"] == sequence].copy()
				start = sequence_set["sstart"].to_string(index = False)
				end = sequence_set["send"].to_string(index = False)
				if blast_re[(blast_re["qstart"] == int(start))].empty:
					if (blast_re["qend"] == int(end)).empty:
						LVL5.append(0)
					else:
						LVL5.append(0)
				else:
					LVL5.append(0.5)
	se = pd.Series(LVL5)
	numts_LVL4['LVL5'] = se.values
	numts_LVL5 = numts_LVL4.copy()
	return(numts_LVL5)

def numts_max(numts_LVL5):
	levels = numts_LVL5["LVL1"].astype(int) + numts_LVL5["LVL2"].astype(int) + numts_LVL5["LVL3"].astype(int) + numts_LVL5["LVL4"].astype(float) + numts_LVL5["LVL5"].astype(float)
	numts_LVL5["LVL_MAX"] = levels
	numts_LVL5['evalue'] = numts_LVL5['evalue'].astype('float64')
	numts_LVL5['pident'] = numts_LVL5['pident'].astype('float64')
	return numts_LVL5
