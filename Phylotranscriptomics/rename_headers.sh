###########################
## bash rename_headers.sh -e [extension]
## Rename taxa headers
## Jessie Pelosi
## Last updated September 29, 2021
###########################

#!/bin/sh 

while getopts e: flag
do
    case "${flag}" in
        e) ext=${OPTARG};;
    esac
done

## 1KP changes --> 
## 		ORJE --> Phymatosorus_grossus_ORJE
##		AFPO --> Athyrium_sp_AFPO 
##		BMJR --> Adiantum_raddianum_BMJR
##		DCDT --> Gaga_arizonica_DCDT 
##		GANB --> Aslophila_spinulosa_GANB 
##		ZQYU --> Phlebodium_pseudoaureum_ZQYU 
##		GSXD --> Myriopteris_rufa_GSXD

for file in *.$ext; do
        sed -i -r 's/scaffold\-([A-Z]*)\-[0-9]*\-([A-Za-z\_\-]*)\.[a-z0-9]*/\2_\1/g' "$file"
        sed -i -r 's/_scaffold[0-9]*\.[a-z0-9]*//g' "$file"
        sed -i -r 's/_C[0-9]*\.[a-z0-9]*//g' "$file"
        sed -i -r 's/Gb_[0-9]*/Ginkgo_biloba/g' "$file"
        sed -i -r 's/AT[0-9]G[0-9]*/Arabidopsis_thaliana/g' "$file"
        sed -i -r 's/AMTR_[A-Za-z0-9]*/Amborella_trichopoda/g' "$file"
        sed -i -r 's/SELMODRAFT_[0-9]*/Selaginella_moellendorffii/g' "$file"
        sed -i -r 's/Pp3[A-Za-z0-9\_]*/Physcomitrella_patens/g' "$file"
        sed -i 's/Azolla_cf_CVEG_caroliniana.p[0-9]/Azolla_caroliniana_CVEG/g' "$file"
        sed -i 's/Rhachidosorus\_sp\.XQ\-2018\_RHSP\_scaffold[0-9]*\_[A-Z0-9]*\.p[0-9]/Rhachidosorus\_sp\.XQ\-2018_RHSP/g' "$file"
        sed -i 's/Osmunda\_sp\_UOMY\-gametophyte\.p[0-9]/Osmunda\_sp\_UOMY\-gametophyte/g' "$file"
        sed -i 's/Hymenophyllum\_sp\.\_XQ\-2018\_HYME\_[A-Za-z0-9\_]*\.p[0-9]/Hymenophyllum\_sp\.\_XQ\-2018\_HYME/g' "$file"
        sed -i 's/Cibotium_glaucum_ORJE/Phymatosorus_grossus_ORJE/g' "$file"
        sed -i 's/Blechnum_spicant-gametophyte_AFPO/Athyrium_sp_AFPO/g' "$file"
        sed -i 's/Adiantum_tenerum_BMJR/Adiantum_raddianum_BMJR/g' "$file"
        sed -i 's/Cheilanthes_arizonica_DCDT/Gaga_arizonica_DCDT/g' "$file"
        sed -i 's/Cyathea_spinulosa_GANB/Alsophila_spinulosa_GANB/g' "$file"
        sed -i 's/Polypodium_plectolens_ZQYU/Phlebodium_pseudoaureum_ZQYU/g' "$file"
        sed -i 's/Myriopteris_eatonii_GSXD/Myriopteris_rufa_GSXD/g' "$file"
        sed -i 's/Pleopeltis_polypodioides-dehydrating_fronds_UJWU/Pleopeltis_polypodioides_UJWU/g' "$file"
        sed -i 's/Pteris_ensigormis_FLTD/Pteris_ensiformis_FLTD/g' "$file"
done
