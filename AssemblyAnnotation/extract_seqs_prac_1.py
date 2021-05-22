'''
extract_seqs_prac_1.py

extract_seqs_prac_1.py [transcrtipts.txt] [transcriptome.fasta] 

Purpose: Extract sequences corresponding to the transcripts file generated using sed during the Trinotate pipeline. 

Written By: Jessie Pelosi
Last Modified: May 22, 2021 

'''


import sys
args = sys.argv 
gene_list = sys.argv[1] 
fasta_name = sys.argv[2]
#for arg in args[2:len(args)]: #use this to search through multiple fasta files sequentially 
#   fasta_names.append(arg)

#from Bio import SeqIO

headers=[] # make a vector of identifers from the input file 
with open(gene_list, "rt") as genes:
	for gene in genes:
		headers.append(gene.strip()) # add each identifer to the vector 

print(headers)

records_to_export = [] # make vector for sequence records 
for header in headers:
	with open(fasta_name, "rt") as fasta:
		for line in fasta:
			if line.startswith(">"):
				flag = False
			if header in line.strip():
				flag = True
			if flag:
				records_to_export.append(line)
'''
lines_to_export = [] # store headers and sequence lines in this vector
for contig in contigs:
    # open each transcriptome file and loop through its contents
    for txm_name in txm_names:
        with open(txm_name, "rt") as txm_file:
            for line in txm_file:
                if line.startswith(">"):
                    flag = False # turns off signal to export sequence once loop hits a new header.
                                 # even the header that matches will trip this if statement, but it'll
                                 # immediately trip the next one too, which turns the signal on.
                if line.strip() == contig: # don't forget to strip whitespace from lines in txm file before comparing
                    flag = True # turns on signal to export once loop hits header that matches
                if flag: # if signal is on, add line to storage (this allows for multi-line sequences)
                    lines_to_export.append(line)
'''              
#write sequences out to file 
with open(fasta_name + "_matching.fasta", "w") as output_file:
    for export in records_to_export:
        output_file.write(export)
