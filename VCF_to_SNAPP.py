#!/usr/bin/python

'''
takes as input a VCF file with one SNP per locus from ipyrad (you have to call one snp per locus somehow)
converts to SNAPP input file, where 0 is homoz allele 1, 1 is hetero, 2 is homoz allele 2
SNAPP requires that SNPs are biallelic
most of my comments assume a dataset with 73 individuals and 7066 loci. Obviously yours would be different.

usage:
VCF_to_SNAPP.py -i input.vcf

written by Kyle O'Connell and Blake O'Connell
oconnellk@si.edu
April 2019
'''

import argparse
import os
import subprocess as sp
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np

#global lists
indlist = [] #list of indiv names from the vcf header
allele_list = [] #large matrix of converted alleles, indiv seq is vertical, so in my data it is 73 ind wide X 7066 loci high
sequence_list = [] # 'flipped' so speak of allele_list, so each indiv seq is horiz , so 7066 loci wide X 73 ind high

#define input variables for input file and out dir
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_vcf", required=True, help="REQUIRED: The full path to the input vcf")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: Specify the output directory")
    return parser.parse_args()

#open input vcf
def file_reader(vcf):
	fh = open(vcf,'r')
	for line in fh:
		#skip extra header rows
		if line.startswith('##'):
			pass
		#grab header with the ind names
		elif line.startswith('#'):
			#grab the header line, split by tab, skip the first 8 cols until you get to the indiv IDs
			for item in line.split('\t')[9:]:
				item = item.strip('\n')
				#append to ind list
				indlist.append(item)
				#create empty list that we will append to later to flip the vertical allele array to horiz array
				sequence_list.append('')
	#close vcf file
	fh.close()

#function to recode ./., 0/0, 0/1, 1/0, 1/1 allele codes in vcf file to 0,1,2,- of SNAPP input
#needs vcf as input
def seq_assign(vcf):
	#open the file again
	fh = open(vcf,'r')
	for line in fh:
		#skip all headers and go to sequence lines
		if line.startswith('#'):
			pass
		else:
			#create empty list as we iterate through each line
			alleles = []
			#homo 1 = 0, homo 2 = 2, het = 1
			#skip first 8 cols and get to indiv alleles
			for allele in line.split('\t')[9:]:
				#grab only allele info, remove location info
				allele = allele.split(':')[0]
				#missing locus
				if allele == './.':
					#replace ./. with - for snapp missing value
					allele = allele.replace('./.','-')
					#append to the alleles list inside this loop which only has the current locus alleles
					#the rest of this section is the same but for each allele value, 0/0, 0/1, 1/0, 1/1
					alleles.append(allele)
				#homoz 1
				elif allele == '0/0':
					allele = allele.replace('0/0','0')
					alleles.append(allele)
				#hetero 
				elif allele == '1/0':
					allele = allele.replace('1/0','1')
					alleles.append(allele)
				#hetero
				elif allele == '0/1':
					allele = allele.replace('0/1','1')
					alleles.append(allele)
				#homoz 2
				elif allele == '1/1':
					allele = allele.replace('1/1','2')
					alleles.append(allele)
			#append this locus line to the large allele array in the global variables at top
			allele_list.append(alleles)
			
#now we convert the large global array of alleles and flip it so that it is horizontal for each individual
def write_sequences():
	#here i is a number, 1-7066, since that is the number of loci I have, in this case organized vertically
	for i in range(len(allele_list)):
		#create a temp list that is just the current line of values, in this case 73 allele values long
		allele_array=allele_list[i]
		#r is 0-72, we want to grab each allele value in that line, and append it to the blank list we made above
		#thus we grab value one which pertains to ind 1, and we write it to the sequence list
		for r in range(len(allele_array)):
			#we have to make these loops so that it can append correctly
			#Blake O'Connell helped with a lot of this part
			if allele_array[r] == '0':
				sequence_list[r] = sequence_list[r] + allele_array[r]
			elif allele_array[r] == '1':
				sequence_list[r] = sequence_list[r] + allele_array[r]
			elif allele_array[r] == '2':
				sequence_list[r] = sequence_list[r] + allele_array[r]
			if allele_array[r] == '-':
				sequence_list[r] = sequence_list[r] + allele_array[r]

#now write a file
#the important part is the loop, where it is matching the line of the indlist with the line of the sequence_list				          
def FileWriter():
    filename2 = "snapp_input_vcf.nex"
    nexusFile = open(filename2, 'w')
    
    nexusFile.write("#NEXUS\n")
    nexusFile.write("Begin data;\n")
    nexusFile.write("\tDimensions ntax=" + str(len(indlist)) + " nchar=" + str(len(allele_list)) + ";\n")
    nexusFile.write("\tFormat datatype=integerdata symbols=\"012\" gap=\"-\";\n")
    nexusFile.write("\tMatrix\n")
    for i in range(len(indlist)):
        nexusFile.write(indlist[i] + "\t" + sequence_list[i] + "\n")
    nexusFile.write("\t;\n")
    nexusFile.write("End;")
    print "output SNAPP input file called {}".format(filename2)
    nexusFile.close()
    

#call all the functions    
def main():
	#define the arguments
	args = get_args()
	file_reader(args.in_vcf)
	seq_assign(args.in_vcf)
	write_sequences()
	FileWriter()
	


if __name__ == '__main__':
    main()




