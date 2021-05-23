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

Primary transcripts for each outgroup were extracted using primary_transcripts.py (from XX) and used in downstream analyses. 
```
python primary_transcripts.py *.cds
python primary_transcripts.py *.pep
```


## References 

Guan R., Y. Zhao, H. Zhang, G. Fan, X. Liu; W. Zhou; C. Shi, J. Wang, W. Liu, X. Liang, Y. Fum, K. Ma, L. Zhao, 
F. Zhang, Z. Lu, S.M. Lee, X. Xu, J. Wang, H. Yang, C. Fu, S. Ge, W. Chen (2019) [Updated genome assembly of Ginkgo biloba.](http://gigadb.org/dataset/100613) GigaScience Database. 
http://dx.doi.org/10.5524/100613
