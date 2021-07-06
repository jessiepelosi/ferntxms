# Repeats

## Repeat Annotation

Repeats were identified and annotated in each transcriptome using [RepeatMasker ver. 4.0.5](https://www.repeatmasker.org/) (Smit et al. 2015) with the viridiplantae repeat database. 
```
RepeatMasker Transcriptome.fasta -pa 6 -species viridiplantae 
```
## Extract Repeat Sequences

The output from RepeatMasker was parsed with the perl pipeline [One Code to Find Them All](http://doua.prabi.fr/software/one-code-to-find-them-all) (Bailly-Bechet et al. 2014)
```
perl buld_dictionary.pl --rm RepeatMasker.out > Txm.LTR
perl one_code_to_find_them_all.pl --rm RepeatMasker.out --ltr Txm.LTR --fasta Txm.fa 
```
Note: Txm.fa needs to be the file that was used as input for RepeatMasker 

The OCTFTA output was then passed to the custom bash script, ThereAndBackAgain.sh. 
```
for line in $(cat subfamily_list.txt); do bash ThereAndBackAgain.sh -r “$line” -s N -p taxon_name; done 
```

You’ll get a lot of files- there are around 7000 TE subfamilies, so there should be one fasta file per subfamily per transcriptome. There will be a bunch of subfamilies where the sequence files are empty- that is okay! That simply means that that particular subfamily was not found in the transcriptome 
We then have to put the sequence files for each subfamily into a single file for all the transcriptomes, then we can align them and make a tree. 

The sequences from all transcriptomes for each TE subfamily were then concatenated using the following: 
```
for line in $(cat subfamilies_list.txt); do cat ./*/*”$line” > txms_“$line”.fasta 
```
This command assumes that the file structure is: current_directory/txm/out_files.fasta


## Repeat Phylogenetics 


## References

Bailly-Bechet M, Haudry A, Lerat E (2014) “One code to find them all”: a perl tool to conveniently parse RepeatMasker output files. Mob DNA 5:13

Smit, A.F.A, R. Hubley, and P. Green. 2013-2015. RepeatMasker Open-4.0 http://www.repeatmasker.org.
