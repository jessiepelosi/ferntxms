
#####################################
# get_orthogroups.sh
# Purpose: to create a file that contains the names of the OGs that contain X number of taxa that are single copy (-s Y)
# or multicopy
# bash get_OGs.sh -l [list_of_taxon_identifiers.txt] -f [OG.fasta] -s [Y/N] -p [0 to 1] 
#
# Jessie Pelosi
# Last modified August 2, 2021
#####################################

#!/bin/bash
while getopts l:f:s:p: flag
do
    case "${flag}" in
        l) taxon_list=${OPTARG};;
        f) OG_file=${OPTARG};;
        s) single_copy=${OPTARG};;
        p) percent=${OPTARG};;
    esac
done

num_taxa=$(wc -l < "$taxon_list")  

if single_copy=Y; then

        for taxon in $(cat $taxon_list);
                 do declare "$taxon"=$(grep -o "$taxon" $OG_file | wc -w) 
                 if [ $(echo "${!taxon}") -gt 1 ]; then declare "$taxon"=para;fi 
                 echo "${!taxon}" >> "$OG_file"_temp.txt
        done

        while read num; do ((sum += num)); done < "$OG_file"_temp.txt;
        proportion=$(echo $sum / $num_taxa | bc -l)
        if (( $(echo "$proportion > $percent" | bc -l) )); then echo "$OG_file" >> SC_OGs_with_"$percent"_proportion.txt;fi


else
        for taxon in $(cat "$taxon_list"); do  
                 declare "$taxon"=$(grep -o "$taxon" $OG_file | wc -w) 
                 if [ $(echo "${!taxon}") -gt 0 ]; then declare "$taxon"=1;fi 
                 echo "${!taxon}" >> "$OG_file"_temp.txt
        done

        while read num; do ((sum += num)); done < "OG_file"_temp.txt;
        proportion=$(echo $sum / $num_taxa | bc -l)
        if (( $(echo "$proportion > $percent" | bc -l) )); then echo "$OG_file" >> MC_OGs_with_"$percent"_proportion.txt;fi
        rm "$OG_file"_temp.txt

fi 
