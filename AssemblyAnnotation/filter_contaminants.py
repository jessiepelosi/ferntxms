'''
filter_contaminants.py 

Usage:

python filter_contaminants.py [input.fasta] [input_blast] > [output_fasta]

arguments: input.fasta = transcriptome assembly
		   input_blast = output from blastn in outfmt6
		   outout_fasta = output transcriptome assembly 

Jessie Pelosi 
Last Modified: May 10, 2021

'''

import sys
from Bio import SeqIO

args = sys.argv 
input_fasta = sys.argv[1]
blast_output = sys.argv[2]

import csv

matches = []
with open(blast_output, "r") as blast:
	rd = csv.reader(blast, delimiter = '\t')
	for row in rd:
		if float(row[3]) > 300:
			if float(row[11]) > 50: 
				matches.append(row[0])


with open (input_fasta + "blast_matches", "w") as blast_hits:
	for s in matches:
		blast_hits.write("%s\n" % s)


header_set = set(line.strip() for line in open(input_fasta + "blast_matches"))

fasta_file = SeqIO.parse(input_fasta, "fasta")

for seq_record in fasta_file:
	try:
		 header_set.remove(seq_record.name)
	except KeyError:
		print(seq_record.format("fasta"))
		continue
