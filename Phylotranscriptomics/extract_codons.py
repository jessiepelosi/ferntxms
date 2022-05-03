'''
Purpose: python code to extract first/second codons third codon positions from fasta file
Usage: python extract_codons.py input_file
Arguments:
    input_file: name of FASTA file from which to extract codons
Last modified May 3, 2022 by Jessie Pelosi
'''
import sys

args = sys.argv
infile_name = sys.argv[1]

from Bio import SeqIO

thirdpos = []
firstsecondpos = []

with open(infile_name, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seqid = '>' + record.description
        seq = str(record.seq)
        cp3 = []
        for start in range(0,len(record.seq),3):
            x = seq[start:start+3]
            if 3 <= len(x):  
                base = x[2] 
                cp3.append(base)
        cp3 = ''.join(cp3)
        posthree = '\n'+ seqid + '\n' + cp3
        thirdpos.append(posthree)
        cp12 = []
        for start in range(0,len(record.seq),3):
            x = seq[start:start+3]
            if 1 <= len(x) or 2 <= len(x):
                bases = x[0:2]
                cp12.append(bases)
        cp12 = ''.join(cp12)
        posonetwo = '\n' + seqid + '\n' + cp12
        firstsecondpos.append(posonetwo)
        

with open(infile_name + ".cp3.fna","w") as output_1:
    for id in thirdpos:
        output_1.write(id)

with open(infile_name + ".cp12.fna","w") as output_2:
    for id in firstsecondpos:
        output_2.write(id)
