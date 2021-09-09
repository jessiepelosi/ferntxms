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
We then used the R package [mixtools ver. 1.2.0](https://cran.r-project.org/web/packages/mixtools/mixtools.pdf) (Benaglia et al. 2009) to fit mixture models to the Ks distributions. The code for this can be found in this directory under `ks_plots.R`. 


## Phylogenomic Approaches 

## Biased Gene Retention 

To examine how and if genes were differentially retained after WGDs, we created lists of duplicates that were +/- 1SD from the mean putative WGD peak(s) in each transcriptome (see `ks_plots.R`). The output files were then further edited to get a list of paralogs that could then be used to extract the corresponding CDS from each transcriptome for further analyses. 

```
for file in *.tsv; do sed -r 's/^[0-9]*\s//g' "$file" | sed 's/\_\_/\n/g' | sed -r 's/^([A-Za-z])/>\1/g' > "$file"_paralog.list.tsv;done
```

We then extracted the corresponding CDS sequences with `extract_cds.py`. 

```
python extract_cds.py sample_paralog.list.tsv sample.fa.cds 
```

We then used the GOGetter Pipeline (E.B. Sessa, M.S. Barker et al., unpub.) to BLAST each transcriptome and corresponding paralogs to the Araport11 <i>Arabidopsis thaliana </i> protein dataset (Berardini et al. 2004). 

```
perl 0_Get_GO_annotations.pl taxon_list.txt 
```

The proption of GO terms in each GO category (as defined by Araport) was visualized in R (see `gene_retention.R`). 


## References 

Zwaenepoel, A., and Van de Peer, Y. wgd - simple command line tools for the analysis of ancient whole genome duplications. Bioinformatics., bty915, https://doi.org/10.1093/bioinformatics/bty915
