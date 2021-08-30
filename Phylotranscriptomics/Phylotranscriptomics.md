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

Primary transcripts for each outgroup were extracted using primary_transcripts.py (from David Emms, accessible [here](https://github.com/davidemms/OrthoFinder/blob/master/tools/primary_transcript.py)) and used in downstream analyses. 
```
python primary_transcripts.py *.cds
python primary_transcripts.py *.pep
```
## 2. Run OrthoFinder 

Place all pep sequences into a single folder named "Proteomes". Run [OrthoFinder v.2.3.11](https://github.com/davidemms/OrthoFinder). 
```
orthofinder -M msa -A mafft -T fasttree -t [threads] -a [threads] -f proteomes/
```
Note that OrthoFinder continously failed due to the large number of transcriptomes being analyzed even with 1Tb of RAM. We instead ran Diamond Blasts independently after the commands were generated using: 
```
orthofinder -f proteomes/ -op
```
We then ran each Diamond Blast in parallel. Once the Blasts were completed, we ran:
```
orthofinder -a 30 -b proteomes/OrthoFinder/Results_XXX/WorkingDirectory/
```

## 3. Filter Orthologs 

We used the custom bash script <i>get_orthogroups.sh </i>to filter for single- and multi-copy orthogroups which contained a certain proportion of the input taxa. The input for this script is a list of unique identifiers (denoted by the `-l` flag), the orthogroup fasta file (denoted by the `-f` flag), whether you want single- or multi-copy orthogroups (either `Y` or `N` for the `-s` flag), and the proption of taxa you wish to filter for (0 to 1, denoted by `-p` flag). 

```
bash get_orthogroups.sh -l [list_of_taxon_identifiers.txt] -f [OG.fasta] -s [Y/N] -p [0 to 1]
```

We first created a set of single-copy orthogroups for which at least 80% of the taxa were present. In the `OrthogroupSequences` directory of the OrthoFinder results, we ran:
```
for file in *.fa; do bash -l unqiue_IDs.txt -f "$file" -s Y -p 0.8;done
```

We also created a set of multi-copy orthogroups for which all taxa were present (100%). In the `OrthogroupSequences` directory of the OrthoFinder results, we ran:
```
for file in *.fa; do bash -l unique_IDs.txt -f "$file" -s N -p 1;done
```

## 4. Extract Corresponding coding sequences 

Since OrthoFinder takes peptide sequences as input, we had to extract the corresponding coding sequences (CDS) for each orthogroup, which we then aligned in Step 5. We used the custom python script <i>extract_CDS.py</i> (modified from Kasey K. Pham). The input_file is the name of the FASTA file to be converted, in this case this is the orthogroup fasta from OrthoFinder. The transcriptome files are the CDS files that are output from trinnotate (for example). 

```
python extract_cds.py input_file trancriptome_1 transcriptome_2 ... transciptome_n
```

This script will only work if the headers in the orthogroup fasta file (and by extension the proteome files) and the transcriptome (CDS) files. Run as loop or array on SLURM. 


## 5. Sequence Alignment + Trimming 

The codon-aware alignment program [MACSE ver. 2.04](https://bioweb.supagro.inra.fr/macse/) (Ranwez et al. 2011) was used to align the extracted CDS sequences for each orthogroup. Run as loop or array on SLURM. 
```
for file in *.fa; do java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq "$file"; done
```

Since we are using so many transcriptomes over deep time the raw alignments from MACSE are a bit dirty. We need to clean these up by removing gappy sites - those that contain less than 50% of the transcriptomes. We can do this with [TrimAl ver. 1.2](http://trimal.cgenomics.org/) (Capella-Gutierrez et al. 2009). 

```
trimal -in OGXXXX_NT.fa -out OGXXX_NT.trim50.fa -gt 0.5
```

In some cases, trimming led to problems where some sequences were composed entirely of gaps. We learned of this issue by running into the log error messages in IQTREE. This can be done in the following way with [qiime ver. 1.9.1](http://qiime.org/) (Caporaso et al. 2010). 

```
for file in $(cat for_IQTREE.txt); do
         grep "WARNING" "$file".log | sed -r 's/WARNING\:\sSequence\s([A-Za-z0-9\_\.\-]*).*/\1/g' > "$file"_seq-to-remove.txt
         filter_fasta.py -f "$file" -o "$file".seqsremoved.fasta -s "$file"_seq-to-remove.txt -n
done
```

The output files from this program are aligned CDS and peptide files. We will use both for gene tree construction (next step). 

## 5. Gene Tree Construction 

Gene trees were constructed under a maximum likelihood framework with [IQTREE2 ver. 2.1.2](http://www.iqtree.org/) (Minh et al. 2020) with ModelFinder (Kalyaanamoorthy et al. 2017) and 1000 ultrafast bootstraps (Hoang et al. 2018). Given the results of Shen et al. (2020) we ran IQTREE2 with two threads for all jobs. Run as loop or array.    

```
iqtree2 -s OGXXX_NT.trim50.fa -m TEST --alrt 1000 -B 1000 -T 2 --redo
```
Run with both peptide and CDS alignments from MACSE. 

## 6. Generate Species Tree 

Prior to running ASTRAL, we must change names of the tips used in the phylogeny so that each transcriptome is only represented once, not per scaffold. 

```
cat *.treefile > MC_loci.txt # Create a single file with every gene tree for the multi-copy dataset
cat *.treefile > SC_loci.txt # Create a single file with every gene tree for the single-copy dataset 
bash rename_trees.sh # Rename the scaffolds to represent their corresponding transcriptome 
```
Then, we can run ASTRAL. 

```
ASTRAL Pro and ASTRAL 
```

We can also compare the multi-species coalescent tree from ASTRAL to a concatenated analysis with the single-copy orthogroup dataset. Orthogroup alignements can be concatenated in Geneious, a partition file created, and run under maximum likelihood with IQTREE 2 as above. 
```
Concat IQTREE2 commands 
```

## 7. Divergence Time Estimation 

```
MARE + LBScore.py 
MCMCTREE
```

## References 

Emms, D.M., and S. Kelly. 2015. Solving fundamental biases in whole genome comparisons drastically improves orthogroup inference accuracy. Genome Biology 16: 157.  

Emms, D.M, and S. Kelly. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology 20: 238. 

Guan R., Y. Zhao, H. Zhang, G. Fan, X. Liu; W. Zhou; C. Shi, J. Wang, W. Liu, X. Liang, Y. Fum, K. Ma, L. Zhao, 
F. Zhang, Z. Lu, S.M. Lee, X. Xu, J. Wang, H. Yang, C. Fu, S. Ge, W. Chen (2019) [Updated genome assembly of Ginkgo biloba.](http://gigadb.org/dataset/100613) GigaScience Database. 
http://dx.doi.org/10.5524/100613

Hoang, D.T., O. Chernomor, A. von Haeseler, B.Q. Minh, and L.S. Vinh. 2018. UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol. 35: 518-522.

Kalyaanamoorthy, S., B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin. 2017. ModelFinder: Fast model selection for accurate phylogenetic estimates. Nat Methods 14: 587-589. 

Meyer, B. and Misof, B  2010 MARE: MAtrix REduction - A tool to select optimized data subsets from supermatrices for phylogenetic inference.. Zentrum für molekulare Biodiversitätsforschung (zmb) amZFMK, Adenauerallee 160, 53113 Bonn, Germany. 

Minh, B.Q., H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534

Puttick, M.N. 2019. MCMCtreeR: functions to prepare MCMCtree analyses and visualize posterior ages on trees. Bioinformatics 35(24): 5321-5322. 

Vincent Ranwez, Sébastien Harispe, Frédéric Delsuc, Emmanuel JP Douzery. MACSE: Multiple Alignment of Coding SEquences accounting for frameshifts and stop codons. PLoS One 2011, 6(9): e22594.

Yang, Z. 2007. PAML 4: a program package for phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24: 1586-1591. 

Zhang, C., C. Scornavacca, E.K. Molloy, and S. Mirabab. 2020. ASTRAL-Pro: Quartet-based species-tree inference despite paralogy. Molecular Biology and Evolution msaa139. 
