# ThereAndBackAgain.sh 
# Usage: bash ThereAndBackAgain.sh -r [repeat_family]  -s [get list of subfamilies? Y/N] -p [output prefix, taxon name] 
# Purpose: Used to extract sequences from output of One Code to Find Them All (CITE); it will aggregate all the sequences for a given TE family/subfamily for a given 
# genome/transcriptome's RepeatMasker output into a single file. For comparative analyses, these files can be concatenated and aligned for tree building. 
# This script can also be used to generate a list of subfamilies given a family of TEs. 
# Examples of Repeat Families: 
#	Gypsy, Copia 
# Examples of subfamilies: 
#       Gypsy-9-I_CR, Gypsy-52_Mad-I
# Last Modified: June 27 2021
# Jessie Pelosi, University of Florida 

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "ThereAndBackAgain.sh is used to further parse the output from "
   echo
   echo "Usage: bash ThereAndBackAgain.sh -r -s -p -h"
   echo "options:"
   echo "r     Repeat family or subfamily."
   echo "s     Get list of subfamilies? [Y/N]."
   echo "p     Output prefix, preferably taxon name."
   echo "h     Print this help message."
   echo
}

###############################################################################
# Main Script                                                                 #
###############################################################################

while getopts r:s:t: flag
do
    case "${flag}" in
        r) repeat_family=${OPTARG};;
	s) subfam=${OPTARG};;
        t) taxon_name=${OPTARG};;
    esac
done

#move all .fasta files into one directory
if [[ -d fastas ]]
then
	echo "fastas/ already exists!"
else
	mkdir fastas
	mv *.fasta fastas/
fi

#move all .csv files into one directory
if [[ -d csvs ]]
then
	echo "csvs/ already exists!"
else
	mkdir csvs
	mv *.csv csvs/ 
fi 

#genreate file with all sequences for specified family of rpeats
repeat_family_files=$(grep -i "$repeat_family" fastas/*.fasta -l)

if [[ -d $repeat_family ]]
then
	echo "$repeat_family already exists!"
else
	mkdir $repeat_family
fi 

for file in $repeat_family_files; do cp "$file" $repeat_family;done

cat $repeat_family/*.fasta > "$repeat_family"_all_seq.temp

grep -i -A 1 "$repeat_family" "$repeat_family"_all_seq.temp > "$taxon_name"_"$repeat_family".fasta 

rm *.temp 

#genreate list of unique subfamilies within the specified family of repeats 
if [ $subfam = Y ]
then 
	for line in $(ls fastas/*.fasta); do grep "$repeat_family" "$line" -h | sed -r 's/.*\|.*\|.*\|([.]*)/\1/g' >> "$repeat_family"_subfamilies.txt;done
	sort "$repeat_family"_subfamilies.txt | uniq > "$repeat_family"_subfamilies.txt
else
	echo "List of subfamilies not reported."
fi 

