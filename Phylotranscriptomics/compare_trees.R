#########################
## compare_trees.R
## Jessie Pelosi
## Last modified December 16, 2021
##
#########################

library(ape)
library(phangorn)
library(pheatmap)

treelist <- read.tree("../OrthoFinder_8.2.21_new/SpeciesTrees/SpeciesTrees_10.1.21/species_tree_list.txt")

distance <- RF.dist(treelist)
normalized <- RF.dist(treelist, normalize = T)

mat <- as.matrix(distance)
normmat <- as.matrix(normalized)

annotcol <- data.frame(Analysis = (c("ASTRAL","ASTRAL","ML","ML","ASTRAL","ASTRAL","ML","ML","ASTRAL","ASTRAL",
                                     "ML","ML", "ASTRAL","ASTRAL", "ASTRAL.GHOST.Unlinked", "ASTRAL.GHOST.Linked", "ASTRAL.CP12", "ASTRAL.CP3")))
annotrow <- data.frame(Dataset = c("SCO60","SCO60","SCO60","SCO60", "SCO75","SCO75","SCO75","SCO75", 
                   "SCO85","SCO85","SCO85","SCO85", "MCO","MCO", "SCO60", "SCO60", "SCO60", "SCO60"))
annotrow$DataType <- c("AA", "NT", "AA", "NT","AA", "NT","AA", "NT","AA", "NT","AA", "NT", "AA","NT", "NT", "NT", "NT", "NT")
row.names(annotcol) <- colnames(normmat)
row.names(annotrow) <- colnames(normmat)

pheatmap(mat, display_numbers = T)
pheatmap(normmat, display_numbers = T, annotation_col = annotcol, annotation_row = annotrow)
