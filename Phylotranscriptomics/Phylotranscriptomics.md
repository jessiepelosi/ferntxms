# Phylotranscriptomics

## 1. Download Outgroups

Outgroups were downloaded from Ensembl Plants 51. We also downloaded the coding sequences and peptide sequences of <i>Ginkgo biloba</i> 
from [Guan et al. 2019](http://gigadb.org/dataset/100613). Outgroups used in this study were: 

|Species                          | Assembly ID|Representative Lineage|
|---------------------------------|------------|----------------------|
|<i>Amborella trichopoda</i>      |AMTR1.0     | Angiosperms          |
|<i>Arabidopsis thaliana</i>      |TAIR10      | Angiosperms          |
|<i>Ginkgo biloba</i>             |NA          | Gymnosperms          |
|<i>Physcomitrium patens</i>      |Phypa V3    | Bryophytes           |
|<i>Selaginella moellendorffii</i>| v1.0       | Lycophytes           |

Primary transcripts for each outgroup were extracted using `primary_transcripts.py` (from David Emms, accessible [here](https://github.com/davidemms/OrthoFinder/blob/master/tools/primary_transcript.py)) and used in downstream analyses. 
```
python primary_transcripts.py *.cds
python primary_transcripts.py *.pep
```
## 2. Run OrthoFinder 

Place all pep sequences into a single folder named "proteomes". Run [OrthoFinder v.2.3.11](https://github.com/davidemms/OrthoFinder). 
```
orthofinder -M msa -A mafft -T fasttree -t [threads] -a [threads] -f proteomes/
```
Note that OrthoFinder failed due to the large number of transcriptomes being analyzed even with 1Tb of RAM. We instead ran Diamond Blasts independently after the commands were generated using: 
```
orthofinder -f proteomes/ -op
```
We then ran each Diamond Blast in parallel using several SLURM arrays. Once the Blasts were completed, we ran:
```
orthofinder -a 30 -b proteomes/OrthoFinder/Results_XXX/WorkingDirectory/
```

## 3. Filter Orthologs 

We used the custom bash script `get_orthogroups.sh` to filter for single- and multi-copy orthogroups which contained a certain proportion of the input taxa. The input for this script is a list of unique identifiers (denoted by the `-l` flag), the orthogroup fasta file (denoted by the `-f` flag), whether you want single- or multi-copy orthogroups (either `Y` or `N` for the `-s` flag), and the proption of taxa you wish to filter for (0 to 1, denoted by `-p` flag). 

```
bash get_orthogroups.sh -l [list_of_taxon_identifiers.txt] -f [OG.fasta] -s [Y/N] -p [0 to 1]
```

We created datasets where 60% (SCO60), 75% (SCO75), and 85% (SCO85) of the taxa were present and single copy. In the `OrthogroupSequences` directory of the OrthoFinder results, we ran:
```
for file in *.fa; do bash -l unqiue_IDs.txt -f "$file" -s Y -p 0.8;done
```

We also created a set of multi-copy orthogroups for which all taxa were present (100%). In the `OrthogroupSequences` directory of the OrthoFinder results, we ran:
```
for file in *.fa; do bash -l unique_IDs.txt -f "$file" -s N -p 1;done
```

## 4. Extract Corresponding coding sequences 

Since OrthoFinder takes peptide sequences as input, we had to extract the corresponding coding sequences (CDS) for each orthogroup, which we then aligned in Step 5. We used the custom python script `extract_CDS.py` (modified from Kasey K. Pham). The input_file is the name of the FASTA file to be converted, in this case this is the orthogroup fasta from OrthoFinder. The transcriptome files are the CDS files that are output from trinnotate (for example). 

```
python extract_cds.py input_file trancriptome_1 transcriptome_2 ... transciptome_n
```

This script will only work if the headers in the orthogroup fasta file (and by extension the proteome files) and the transcriptome (CDS) files. Run as loop or array on SLURM. 

## 5. Sequence Alignment + Trimming 

Prior to aligning sequences for the SCO datasets, we ran `rename_headers.sh`. 
```
bash rename_headers.sh -e cds
```
This script renames the headers of the CDS files to reflect just the taxon name, not individual contigs or scaffolds. As there were always some taxa that contained multiple copies in these mostly single copy datasets, we used [SeqKit ver. 0.10.2](https://bioinf.shenwei.me/seqkit/) (Shen et al. 2016) to remove these duplicates. 
```
for file in *.cds; do 
	seqkit rmdup -n -i "$file" -o "$file"_nodup.fasta
done
```
The codon-aware alignment program [MACSE ver. 2.04](https://bioweb.supagro.inra.fr/macse/) (Ranwez et al. 2011) was used to align the extracted CDS sequences for each orthogroup. Run as loop or array on SLURM. Some alignments required additional memory, up to 20Gb. 
```
for file in *.fa; do java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq "$file"; done
```

There are some characters that MACSE introduces that we need to change before continuing. 
```
for file in *_NT.fasta; do sed -i 's/\!/\-/g' "$file";done
for file in *_AA.fasta; do sed -i 's/\!/\-/g' "$file" | sed -i 's/\*/\-/g';done 
```

Prior to trimming the sequences, we should run `extract_codons.py` to pull out the first/second and third codon positions. We must do this prior to trimming the alignments because TrimAl is not codon-aware and can break potentially break-up codons. `extract_codons.py` takes an codon-aware alignment and outputs two `.fna` files; one for the first and second codon positions and one for the third codon positions. These files can then be passed to TrimAl to remove gappy site prior to building a tree. 

```
for file in *.fasta; do python extract_codons.py "$file";done
```

Since we are using so many transcriptomes over deep time the raw alignments from MACSE are a bit dirty. We need to clean these up by removing gappy sites - those that contain less than 50% of the transcriptomes. We can do this with [TrimAl ver. 1.2](http://trimal.cgenomics.org/) (Capella-Gutierrez et al. 2009). 

```
trimal -in OGXXXX_NT.fa -out OGXXX_NT.trim50.fa -gt 0.5
```

In some cases, trimming led to problems where some sequences were composed entirely of gaps. We learned of this issue by running into the log error messages in IQTREE. This can be done in the following way with Seqkit. 

```
for file in $(cat for_IQTREE.txt); do
         grep "WARNING" "$file".log | sed -r 's/WARNING\:\sSequence\s([A-Za-z0-9\_\.\-]*).*/\1/g' > "$file"_seq-to-remove.txt
         do seqkit grep "$file" -f "$file"_seq-to-remove.txt -v -o "$file".seqsremoved.fasta
done
```

The output files from this program are aligned CDS and peptide files. We will use both for gene tree construction (next step). 

## 6. Gene Tree Construction 

Gene trees were constructed under a maximum likelihood framework with [IQTREE2 ver. 2.1.2](http://www.iqtree.org/) (Minh et al. 2020) with ModelFinder (Kalyaanamoorthy et al. 2017) and 1000 ultrafast bootstraps (Hoang et al. 2018). Given the results of Shen et al. (2020) we ran IQTREE2 with two threads for all jobs. Run as loop or array.    

```
iqtree2 -s OGXXX_NT.trim50.fa -m TEST --alrt 1000 -B 1000 -T 2 --redo
```
Run with both peptide and CDS alignments from MACSE. 

We also ran linked and unliked GHOST models [Crotty et al. 2020](https://academic.oup.com/sysbio/article/69/2/249/5541793) with k=2-6 for the SCO60 NT dataset using the best model identified in the original run: 
```
# Linked
for i in {2..6}; do iqtree -s OGXXX_NT.trim50.fa --alrt 1000 -B 1000 -T 2 --redo -pre OGXXX_NT.trim50.fa.${i} -m [model]+H${i};done
# Unlinked 
for i in {2..6}; do iqtree -s OGXXX_NT.trim50.fa --alrt 1000 -B 1000 -T 2 --redo -pre OGXXX_NT.trim50.fa.${i} -m [model]*H${i};done
```

## 7. Generate Species Tree 

Prior to running ASTRAL, we must change names of the tips used in the phylogeny so that each transcriptome is represented by only its name, not name and scaffold/contig ID. 

```
cat *.treefile > MC_loci.txt # Create a single file with every gene tree for the multi-copy dataset
cat *.treefile > SC_loci.txt # Create a single file with every gene tree for the single-copy dataset 
bash rename_headers.sh -e tre # Rename taxa to represent their corresponding transcriptome 
```
Then, we can run ASTRAL, which uses the mutli-species coalenscent to generate species trees. [ASTRAL ver. 5.7.7](https://github.com/smirarab/ASTRAL) (Zhang et al. 2018) was used to infer species trees with the SCO datasets, and [ASTRAL-Pro ver. 1.1.3](https://github.com/chaoszhang/A-pro) (Zhang et al. 2020) was used to infer species trees for the MCO dataset.

```
#General ASTRAL command (using wrapper script provided by UFRC)
astral -i [gene_tree.tre] -o [species_tree.tre]

#General ASTRAL-Pro command 
java -Djava.library.path=./lib -jar astral.1.1.3.jar -i [gene_trees.tre] -o [species_tree.tre] 
```

We can also compare the multi-species coalescent tree from ASTRAL to a concatenated analysis with the single-copy orthogroup dataset. Orthogroup alignements can be concatenated in Geneious, a partition file created, and run under maximum likelihood with IQTREE2 as above. The partition files are provided in the Data folder. 

For SCO85:
```
iqtree2 -s SCO85_FNA_Concat.fasta --alrt 1000 -B 1000 -p SCO85_FNA_partition.txt --redo -T 2 
iqtree2 -s SCO85_FAA_Concat.fasta --alrt 1000 -B 1000 -p SCO85_FAA_partition.txt --redo -T 2 
```

For SCO75 (Needed to run SCO75 AA with additional CPUs, needed more RAM): 

```
iqtree2 -s SCO75_FNA_Concat.fasta --alrt 1000 -B 1000 -T 2 --redo -p SCO75_FNA_partition.txt
iqtree2 -s SCO75_FAA_Concat.fasta --alrt 1000 -B 1000 -T 10 --redo -p SCO75_FAA_partition.txt
```

For SCO60 (Needed to run with additional CPUs, needed more RAM): 
```
iqtree2 -s SCO85_FNA_Concat.fasta --alrt 1000 -B 1000 -p SCO85_FNA_partition.txt --redo -T 10 
iqtree2 -s SCO85_FAA_Concat.fasta --alrt 1000 -B 1000 -p SCO85_FAA_partition.txt --redo -T 10 
```

Discordance in the data was visualized with [DiscoVista](https://github.com/esayyari/DiscoVista) (Sayyari et al. 2018). Before running DiscoVista, change - to _ 
```
sed -i 's/\-/\_/g' [file]
```
We then need to root the species tree(s) with [newick utils ver. 1.6](https://github.com/tjunier/newick_utils) (Junier and Zdobnov 2010). 
```
nw_reroot [astral tree] Physcomitrella_patens > species_tree_estimated.tree
```
Generate the clade definitions file. I generated two separate files, one for the family-level classification and one for the order-level classification. 
```
python generate_clade-defs.py [annotation file] [outputfile] [Other clades file]
```

Run the discordance analysis on species trees. To use this with bootstraps from IQTREE, need to remove sh-aLRT values. 
```
./discoVista.py -m 0 -c clades-def.txt -p $path -t 95 -o $path/results
```

## 8. Divergence Time Estimation 

We used [SortaDate](https://github.com/FePhyFoFum/SortaDate)(Smith et al. 2018) to generate a list of candidate loci for our divergence time estimation analyses. We used the SCO60 dataset as our input.  

```
# 1. Get root-to-tip variance 
python SortaDate/src/get_var_length.py --flend .treefile --outf var_length.out --outg Physcomitrella_patens .

# 2. Get bipartition support 
python SortaDate/src/get_bp_genetrees.py --flend .treefile --outf bp_support.out . [species_tree.tre]

# 3. Combine the results from these two runs
python SortaDate/src/combine_results.py var_length.out bp_support.out --outf combined.out 

# 4. Sort and get the list of the good genes
python SoratDate/src/get_good_genes.py combined.out --max 3 --order 3,1,2 --outf good_genes.out 
```

We first estimated divergence times with TreePL. Since the ML species tree topology differs from the ASTRAL topology and branch lengths must be in number of substitutions per site, we first generated a maximum likelihood tree with a constrained topology based on the SCO60 FNA ASTRAL topology with IQTREE2:  
```
iqtree2 -s 99Loci_Concat.fasta -m TEST -g 0_SpeciesTree.tre -B 1000 --alrt 1000 --prefix test.1 -p partitions.txt --redo -T 2
```
We then rooted the tree as above and ran TreePL with 15 fossil callibrations. 
```
treePL treePL_15fossils.config
```
We then incorporated phylogenetic uncertainty into the dating analysis following https://github.com/sunray1/treepl. 

We first generated 100 ML BS replicate alignments from the 99 locus alignment. 
```
raxmlHPC-PTHREADS-SSE3 -f j -m GTRGAMMA -b $RANDOM -# 100 -s ../99Loci_Concat.fasta -q ../partitions_4.11.22.txt -n ferns.bs4treepl -T $SLURM_CPUS_PER_TASK
```
Then we generated trees from each alignment with a constrained topology. 
```
for b in {0..99}; do
raxmlHPC-PTHREADS-SSE3 -f e -t ../0_SpeciesTree.tre -m GTRGAMMA -s ../99Loci_Concat.fasta.BS${b} -q ../partitions_4.11.22.txt.BS${b} -n tree_BS${b} -T $SLURM_CPUS_PER_TASK; done
```
Trees were then re-rooted using newick utils: 
```
for file in *.tre; do nw_reroot "$file" Physcomitrella_patens > "$file".reroot;done
```
And dated with TreePL using the same fifteen fossil callibrations. We used the python scripts provided in https://github.com/sunray1/treepl to perform the cross validation and dating analysis. Dated trees were used to annotate the species tree in treeAnnotator in BEAST. 

## Supplemental Analysis 

In a supplemental analysis we downloaded sequence data from three additional taxa. In order to identify orthologs that we could use in generating a new tree, we used BLASTp with each proteome as the query and each orthogroup as the subject. 
```
module load ncbi_blast
for file in */*.fa; do blastp -db "$file" -query [txm_assembly].fa.transdecoder.pep -outfmt 6 -out "$file".[txm] -evalue 0.00001 -max_target_seqs 10;done
# Get best hit based on bit-score 
for file in */*.fa.[txm]; do cat "$file" | sort -rnk 12 | head -n 1 | cut -d$'\t' -f1 > "$file".besthit;done
for dir in */; do cd "$dir"; cat *.besthit > results.fasta; cd ../;done
```
We then used `extract_cds.py` to get the corresponding coding sequence for each orthogroup and used `cat` to generate one file with the new coding sequences from the three additional taxa and the original taxa for each orthogroup. We then aligned and built a species tree as above. 

The we dated the tree using the same methods described above with treePL and python scripts provided in https://github.com/sunray1/treepl for bootstrapping our dating analysis. 
