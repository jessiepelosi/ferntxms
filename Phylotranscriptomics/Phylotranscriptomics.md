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

Place all pep sequences into a single folder named "Proteomes". Run [OrthoFinder v.2.5.2](https://github.com/davidemms/OrthoFinder). 
```
orthofinder -M msa -A mafft -T fasttree -t [threads] -a [threads] -f proteomes/
```

## 3. Filter Orthologs 
```
bash script? 
```

## 4. Sequence Alignment 

```
java -Xms4000m -jar maxse_v2.04.jar -prog alignSequences -seq OGXXXXX.fa 
```

## 5. Gene Tree Construction 

```
iqtree2 -s OGXXXX.NT -m TEST --alrt 1000 -B 1000 -T 16 --redo
```

## 6. Generate Species Tree 

```
ASTRAL Pro and ASTRAL 
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

Yang, Z. 2007. PAML 4: a program package for phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24: 1586-1591. 

Zhang, C., C. Scornavacca, E.K. Molloy, and S. Mirabab. 2020. ASTRAL-Pro: Quartet-based species-tree inference despite paralogy. Molecular Biology and Evolution msaa139. 
