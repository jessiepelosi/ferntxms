# Whole Genome Duplication Analyses

## Ks Plots

Compute Ks plots (paralog age distributions) with [wgd](https://github.com/arzwa/wgd) (Zwaenepoel and Van de Peer 2019). 

Run all-vs.-all Blastp analysis and cluster using MCL with the default parameters (e-value cut-off of 10^-10, inflation factor of 2.0). Make sure that the .cds file is in the working directory!
```
wgd mcl -s sample.fa --cds --mcl
```
Output is a directory named wgd_blast with a file named `sample.fa.blast.tsv.mcl`. Move this file to the directory you’re working with for easy access
```
mv wgd_blast/sample.fa.blast.tsv.mcl /dir_path
```
Compute Ks distributions. This will output a .tsv file and a .svg file. You can check the Ks distribution via the .tsv file and the generated histograms via the .svg file.
```
wgd ksd sample.fa.blast.tsv.mcl sample.fa
```
We then used the R package [mixtools ver. 1.2.0](https://cran.r-project.org/web/packages/mixtools/mixtools.pdf) (Benaglia et al. 2009) to fit mixture models to the Ks distributions. The code for this can be found in this directory under "ks_plots.R". 


## Phylogenomic Approaches 

## References 

Zwaenepoel, A., and Van de Peer, Y. wgd - simple command line tools for the analysis of ancient whole genome duplications. Bioinformatics., bty915, https://doi.org/10.1093/bioinformatics/bty915
