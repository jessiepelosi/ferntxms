'''Purpose: python code to extract coding sequences (cds) from transdecoder output given a fasta file
with the same headers (in particular from OrthoFinder amino acid alignments)
Usage: python extract_cds.py input_file trancriptome_1 transcriptome_2 ... transciptome_n
Arguments:
    input_file: the name of the FASTA file to be converted. Input file is the fasta from OrthoFinder in this case. 
	    You can give a full address here so that
                you don't need to hardcode directories into the script.
    transcriptome_1/2/3: the name of the transcriptome file(s) to search through. Only the first is required
                         to be provided
    ## NOTE: this script is modified from Kasey Pham 
    ## Last modified May 14, 2021
'''

import sys

args = sys.argv # get vector of all arguments passed to command line
infile_name = sys.argv[1] # skips the first argument because that is the script name
# save filename of each transcriptome to scan
txm_names = []
for arg in args[2:len(args)]:
    txm_names.append(arg)

headers = []
with open(infile_name, "r") as infile:
    for line in infile:
        if line.startswith(">"):
            headers.append(line)
    output = ''.join(headers)
    with open(infile_name + ".headers","w") as outfile:
        outfile.write(output); outfile.close()

# save headers in a vector
contigs = []
with open(infile_name + ".headers","rt") as header_file:
    for line in header_file:
        contigs.append(line.strip()) # the strip removes whitespace just in case; you avoid a really 
                                     # frustrating error if there's an extra space or two at the end
                                     # of a header line on one file but not the other.

# loop through all headers and search each transcriptome for each
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

#  export to output file
with open(infile_name + ".cds", "w") as output_file:
    for to_export in lines_to_export:
        output_file.write(to_export)
