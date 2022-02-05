#########################
## MAPS.R
## Jessie Pelosi
## Last Modified Feb 5 2022
##
#########################

setwd("C:/Users/Owner/Dropbox/(Insert cool lab name here)/Jessie/Fern_transcriptomes/fern_txms/")

library(ggplot2)
library(dplyr)
library(ape)
library(ade4)
library(phytools)
library(geiger)
library(WGDgc)
library(reshape2)

#Generate function to calculate geometric mean 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Read in ultrametric tree 
time_tree <- read.tree("../OrthoFinder_8.2.21_new/DivTimeEst/99Loci/Phy/TreePL/treepl_15TEST.tre")

#### EQUI ####
# Trim tree to only include taxa in MAPS analysis 
taxa <- c("Physcomitrella_patens", "Amborella_trichopoda", "Selaginella_moellendorffii",
          "Equisetum_diffusum_EQDI", "Equisetum_hymale_JVSZ",
          "Microsorum_scolopendria_MISC")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "trimmed_tree_toedit.tre")

EQUI_OrthoFinder <-  read.delim("EQUI.Orthogroups.GeneCount")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(EQUI_OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                                'Amborella.trichopoda', 'Microsorum.scolopendria.MISC',
                                'Equisetum.diffusum.EQDI','Equisetum.hymale.JVSZ') %>% 
    filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
             Amborella.trichopoda < 100 & Microsorum.scolopendria.MISC < 100 & 
             Equisetum.diffusum.EQDI < 100 & Equisetum.hymale.JVSZ < 100) %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Microsorum.scolopendria.MISC >=1 | 
          Equisetum.diffusum.EQDI >=1 | Equisetum.hymale.JVSZ >=1 |
          Selaginella.moellendorffii >=1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Microsorum.scolopendria.MISC,
                                Equisetum.diffusum.EQDI, Equisetum.hymale.JVSZ))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.391743 

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("EQUI_test.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)

MLE <- MLEGeneCount(tree, geneCountData = tmp, geomMean = 1.508603, conditioning = "oneInBothClades")

for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.508603,
               conditioning = "oneInBothClades", fixedRetentionRates = T, startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}
birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.0011697
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.00114765 

#### OPHIO.1 #####

taxa <- c("Physcomitrella_patens", "Selaginella_moellendorffii", "Amborella_trichopoda", 
          "Sceptridium_japonicum_SCJA", "Equisetum_arvense_EQAR", "Sceptridium_dissectum_EEAQ", 
          "Ophioderma_pendula_OPPE")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "OPHIO1_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OPHIO.1_OrthoFinder <- read.delim("OPHIO.1.Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OPHIO.1_OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Sceptridium.japonicum.SCJA',
                           'Equisetum.arvense.EQAR','Sceptridium.dissectum.EEAQ', 'Ophioderma.pendula.OPPE') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Sceptridium.japonicum.SCJA < 100 & 
           Equisetum.arvense.EQAR < 100 & Sceptridium.dissectum.EEAQ < 100 &
           Ophioderma.pendula.OPPE < 100)  %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Selaginella.moellendorffii >=1 | 
           Sceptridium.japonicum.SCJA >=1 | Equisetum.arvense.EQAR >=1 |
           Sceptridium.dissectum.EEAQ >=1 | Ophioderma.pendula.OPPE >=1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Sceptridium.japonicum.SCJA,
                        Equisetum.arvense.EQAR, Sceptridium.dissectum.EEAQ, Ophioderma.pendula.OPPE))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.459968

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("OPHIO1_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean =1.459968,
                      conditioning = "oneInBothClades", fixedRetentionRates = T,startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001419779
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001425053

#### OPHIO.2 ####
taxa <- c("Physcomitrella_patens", "Selaginella_moellendorffii", "Amborella_trichopoda", 
          "Sceptridium_japonicum_SCJA", "Equisetum_arvense_EQAR", "Ophioglossum_thermale_OPTH", 
          "Ophioderma_pendula_OPPE", "Ophioglossum_vulgatum_OPVU")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "OPHIO2_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OPHIO.2_OrthoFinder <- read.delim("OPHIO.2.Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OPHIO.2_OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Sceptridium.japonicum.SCJA',
                           'Equisetum.arvense.EQAR', 'Ophioderma.pendula.OPPE','Ophioglossum.thermale.OPTH',
                           'Ophioglossum.vulgatum.OPVU') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Sceptridium.japonicum.SCJA < 100 & 
           Equisetum.arvense.EQAR < 100 & Ophioglossum.thermale.OPTH < 100 &
           Ophioglossum.vulgatum.OPVU < 100)  %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Selaginella.moellendorffii >=1 | 
           Sceptridium.japonicum.SCJA >=1 | Equisetum.arvense.EQAR >=1 |
           Ophioglossum.thermale.OPTH >=1 | Ophioglossum.vulgatum.OPVU >=1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Sceptridium.japonicum.SCJA,
                        Equisetum.arvense.EQAR, Ophioglossum.thermale.OPTH, Ophioderma.pendula.OPPE, Ophioglossum.vulgatum.OPVU))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.446376

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("OPHIO2_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean =1.459968,
                      conditioning = "oneInBothClades", fixedRetentionRates = T,startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001334115
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001378472

#### PSIL ####
taxa <- c("Selaginella_moellendorffii", "Amborella_trichopoda", 
          "Tmesipteris_parva_ALVQ", "Equisetum_arvense_EQAR", "Ophioglossum_thermale_OPTH", 
          "Psilotum_nudum_PSNU", "Tmesipteris_tannensis_TMTA")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "PSIL_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("PSIL_run2_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder, 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Tmesipteris.parva.ALVQ',
                           'Equisetum.arvense.EQAR', 'Ophioglossum.thermale.OPTH','Psilotum.nudum.PSNU',
                           'Tmesipteris.tannensis.TMTA') %>% 
  filter(Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Tmesipteris.parva.ALVQ < 100 & 
           Equisetum.arvense.EQAR < 100 & Ophioglossum.thermale.OPTH < 100 &
           Psilotum.nudum.PSNU < 100 & Tmesipteris.tannensis.TMTA < 100)  %>% 
  filter(Selaginella.moellendorffii >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | 
           Tmesipteris.parva.ALVQ >=1 | Equisetum.arvense.EQAR >=1 |
           Ophioglossum.thermale.OPTH >=1 | Psilotum.nudum.PSNU >=1 | Tmesipteris.tannensis.TMTA >=1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Selaginella.moellendorffii, Amborella.trichopoda, Tmesipteris.parva.ALVQ,
                        Equisetum.arvense.EQAR, Ophioglossum.thermale.OPTH, Psilotum.nudum.PSNU, Tmesipteris.tannensis.TMTA))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.478997

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("PSIL_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean =1.478997,
                      conditioning = "oneInBothClades", fixedRetentionRates = T,startingQ = c(0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # RUN 1: 0.001299071, RUN 2: 0.001536341
death_df <- melt(as.data.frame(death))
mean(death_df$value) # RUN 1: 0.001313255, RUN 2: 0.001437485

#### MARA ####
taxa <- c("Amborella_trichopoda", 
          "Angiopteris_fokiensis_ANFK", "Equisetum_arvense_EQAR", "Psilotum_nudum_PSNU", 
          "Christensenia_aescuilfolira_CHAE", "Ptisana_pellucida_PTPE", "Danaea_nodosa_DANO")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "MARA_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("MARA_run2_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Amborella.trichopoda', 'Angiopteris.fokiensis.ANFK',
                           'Equisetum.arvense.EQAR', 'Psilotum.nudum.PSNU','Christensenia.aescuilfolira.CHAE',
                           'Ptisana.pellucida.PTPE', 'Danaea.nodosa.DANO') %>% 
  filter(Amborella.trichopoda < 100 & Angiopteris.fokiensis.ANFK < 100 & 
           Equisetum.arvense.EQAR < 100 & Psilotum.nudum.PSNU < 100 &
           Ptisana.pellucida.PTPE < 100 & Christensenia.aescuilfolira.CHAE < 100 & Danaea.nodosa.DANO < 100)  %>% 
  filter(Amborella.trichopoda >= 1) %>% 
  filter(Angiopteris.fokiensis.ANFK >=1 | Equisetum.arvense.EQAR >=1 |
           Psilotum.nudum.PSNU >=1 | Ptisana.pellucida.PTPE >=1 | Christensenia.aescuilfolira.CHAE >=1 |
           Danaea.nodosa.DANO >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Amborella.trichopoda, Angiopteris.fokiensis.ANFK,
                        Equisetum.arvense.EQAR, Psilotum.nudum.PSNU, Ptisana.pellucida.PTPE, 
                        Christensenia.aescuilfolira.CHAE, Danaea.nodosa.DANO))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.374465

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("MARA_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.374465,
                      conditioning = "oneInBothClades", fixedRetentionRates = T,startingQ = c(0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001288237
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001146853

#### OSMU ####
taxa <- c("Physcomitrella_patens", "Selaginella_moellendorffii", "Amborella_trichopoda", 
          "Plenasium_banksiifolium_PLBA", "Equisetum_arvense_EQAR", "Microsorum_scolopendria_MISC", 
          "Osmunda_japonica_OSJA", "Osmunda_javanica_VIBO")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "OSMU_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("OSMU.Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Equisetum.arvense.EQAR', 
                           'Microsorum.scolopendria.MISC','Plenasium.banksiifolium.PLBA',
                           'Osmunda.japonica.OSJA', 'Osmunda.javanica.VIBO') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Equisetum.arvense.EQAR < 100 & 
           Microsorum.scolopendria.MISC < 100 & Plenasium.banksiifolium.PLBA < 100 &
           Osmunda.japonica.OSJA < 100 & Osmunda.javanica.VIBO < 100)  %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Selaginella.moellendorffii >=1 | 
           Equisetum.arvense.EQAR >=1 | Microsorum.scolopendria.MISC >=1 |
           Plenasium.banksiifolium.PLBA >=1 | Osmunda.japonica.OSJA >=1 | Osmunda.javanica.VIBO >=1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Equisetum.arvense.EQAR,
                        Microsorum.scolopendria.MISC, Microsorum.scolopendria.MISC, Plenasium.banksiifolium.PLBA, 
                        Osmunda.japonica.OSJA, Osmunda.javanica.VIBO))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.421489

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("OSMU_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.421489,
                      conditioning = "oneInBothClades", fixedRetentionRates = T,startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001255609
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001095926

#### GELI ####
# Trim tree to only include taxa in MAPS analysis 
taxa <- c("Physcomitrella_patens", "Amborella_trichopoda", "Selaginella_moellendorffii",
          "Lygodium_japonicum_LYJA", "Hymenophyllum_sp._XQ_2018_HYME",
          "Sticherus_truncatus_STTR", "Diplopterygium_glucum_DIGL", "Dicranopteris_pedata_DIPE")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "GELI_trimmed_tree_toedit.tre")

OrthoFinder <-  read.delim("GLEI_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Lygodium.japonicum.LYJA',
                           'Hymenophyllum.sp.HYME','Sticherus.truncatus.STTR', 
                           'Diplopterygium.glucum.DIGL', 'Dicranopteris.pedata.DIPE') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Lygodium.japonicum.LYJA < 100 & 
           Hymenophyllum.sp.HYME < 100 & Sticherus.truncatus.STTR < 100 &
           Diplopterygium.glucum.DIGL < 100 & Dicranopteris.pedata.DIPE < 100) %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Lygodium.japonicum.LYJA >=1 | 
           Hymenophyllum.sp.HYME >=1 | Sticherus.truncatus.STTR >=1 |
           Diplopterygium.glucum.DIGL >=1 | Dicranopteris.pedata.DIPE >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Lygodium.japonicum.LYJA,
                        Hymenophyllum.sp.HYME, Sticherus.truncatus.STTR, Diplopterygium.glucum.DIGL,Dicranopteris.pedata.DIPE))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.461235 

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("GELI_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)

for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.508603,
                      conditioning = "oneInBothClades", fixedRetentionRates = T, startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}
birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) #  0.001087586
death_df <- melt(as.data.frame(death))
mean(death_df$value) #  0.0007937658

#### DIPT ####
# Trim tree to only include taxa in MAPS analysis 
taxa <- c("Physcomitrella_patens", "Amborella_trichopoda", "Selaginella_moellendorffii",
          "Lygodium_japonicum_LYJA", "Cheiropleuria_integrifolia_CHIN",
          "Sticherus_truncatus_STTR", "Dipteris_lobbiana_DILO", "Dipteris_conjugata_DICO")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "DIPT_trimmed_tree_toedit.tre")

OrthoFinder <-  read.delim("DIPT_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Lygodium.japonicum.LYJA',
                           'Cheiropleuria.integrifolia.CHIN','Sticherus.truncatus.STTR', 
                           'Dipteris.lobbiana.DILO', 'Dipteris.conjugata.DICO') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Lygodium.japonicum.LYJA < 100 & 
           Cheiropleuria.integrifolia.CHIN < 100 & Sticherus.truncatus.STTR < 100 &
           Dipteris.lobbiana.DILO < 100 & Dipteris.conjugata.DICO < 100) %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Selaginella.moellendorffii >=1 | 
           Lygodium.japonicum.LYJA >=1 | Sticherus.truncatus.STTR >=1 |
           Cheiropleuria.integrifolia.CHIN >=1 | Dipteris.lobbiana.DILO >= 1 |
           Dipteris.conjugata.DICO >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Lygodium.japonicum.LYJA,
                        Cheiropleuria.integrifolia.CHIN, Sticherus.truncatus.STTR, Dipteris.lobbiana.DILO, Dipteris.conjugata.DICO))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.475202

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("DIPT_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)

for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.475202,
                      conditioning = "oneInBothClades", fixedRetentionRates = T, startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}
birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001344581 
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001303272

#### HYMN ####
# Trim tree to only include taxa in MAPS analysis 
taxa <- c("Physcomitrella_patens", "Amborella_trichopoda", "Selaginella_moellendorffii",
          "Lygodium_japonicum_LYJA", "Sticherus_truncatus_STTR", 
          "Cephalomanes_javanicum_CEJA", "Callistopteris_apifolia_CAAP", "Crepidomanes_minutum_CRMI", 
          "Hymenophyllum_sp._XQ_2018_HYME")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "HYMN_trimmed_tree_toedit.tre")

OrthoFinder <-  read.delim("HYMN_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder, 'Physcomitrella.patens', 'Selaginella.moellendorffii',
                           'Amborella.trichopoda', 'Lygodium.japonicum.LYJA',
                           'Sticherus.truncatus.STTR', 'Cephalomanes.javanicum.CEJA',
                           'Hymenophyllum.sp.HYME', 'Callistopteris.apifolia.CAAP', 'Crepidomanes.minutum.CRMI') %>% 
  filter(Physcomitrella.patens < 100 & Selaginella.moellendorffii < 100 &
           Amborella.trichopoda < 100 & Lygodium.japonicum.LYJA < 100 & 
           Sticherus.truncatus.STTR < 100 & Cephalomanes.javanicum.CEJA < 100 &
           Hymenophyllum.sp.HYME < 100 & Callistopteris.apifolia.CAAP < 100 & Crepidomanes.minutum.CRMI < 100) %>% 
  filter(Physcomitrella.patens >= 1) %>% 
  filter(Amborella.trichopoda >= 1 | Selaginella.moellendorffii >=1 | 
           Lygodium.japonicum.LYJA >=1 | Sticherus.truncatus.STTR >=1 |
           Cephalomanes.javanicum.CEJA >=1 | Hymenophyllum.sp.HYME >= 1 |
           Callistopteris.apifolia.CAAP >= 1 | Crepidomanes.minutum.CRMI >= 1) 

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Physcomitrella.patens, Selaginella.moellendorffii, Amborella.trichopoda, Lygodium.japonicum.LYJA,
                        Sticherus.truncatus.STTR, Cephalomanes.javanicum.CEJA, Hymenophyllum.sp.HYME, 
                        Callistopteris.apifolia.CAAP, Crepidomanes.minutum.CRMI))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.38152

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("HYMN_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)

for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.38152,
                      conditioning = "oneInBothClades", fixedRetentionRates = T, startingQ = c(0.2, 0.2, 0.2))
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}
birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001192022
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.0008802873
#### SALV ####
taxa <- c("Azolla_pinnata_AZPN", "Azolla_caroliniana_CVEG", "Salvinia_natans_SANA", "Marsilea_quadrifolia_MAQU", 
          "Lygodium_flexuosum_LYFL", "Sticherus_truncatus_STTR", "Osmunda_japonica_OSJA")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "SALV_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("SALV_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Azolla.pinnata.AZPN', 'Azolla.caroliniana.CVEG',
                           'Salvinia.natans.SANA', 'Marsilea.quadrifolia.MAQU','Lygodium.flexuosum.LYFL',
                           'Sticherus.truncatus.STTR', 'Osmunda.japonica.OSJA') %>% 
  filter(Azolla.pinnata.AZPN < 100 & Azolla.caroliniana.CVEG < 100 & 
           Salvinia.natans.SANA < 100 & Marsilea.quadrifolia.MAQU < 100 &
           Lygodium.flexuosum.LYFL < 100 & Sticherus.truncatus.STTR < 100 & Osmunda.japonica.OSJA < 100)  %>% 
  filter(Osmunda.japonica.OSJA >= 1) %>% 
  filter(Azolla.pinnata.AZPN >=1 | Azolla.caroliniana.CVEG >=1 |
           Salvinia.natans.SANA >=1 | Marsilea.quadrifolia.MAQU >=1 | Lygodium.flexuosum.LYFL >=1 |
           Sticherus.truncatus.STTR >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Azolla.pinnata.AZPN, Azolla.caroliniana.CVEG,
                        Salvinia.natans.SANA, Marsilea.quadrifolia.MAQU, Lygodium.flexuosum.LYFL, 
                        Sticherus.truncatus.STTR, Osmunda.japonica.OSJA))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.253535

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("SALV_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.253535,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001486009
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001307382
#### CYATH1 ####
taxa <- c("Azolla_pinnata_AZPN", "Alsophila_spinulosa_ALSI", "Dicksonia_antarctica_DIAN", "Thyrsopteris_elegans_EWXK", 
          "Lygodium_flexuosum_LYFL", "Sticherus_truncatus_STTR", "Osmunda_japonica_OSJA")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "CYATH1_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("CYATH1_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Azolla.pinnata.AZPN', 'Alsophila.spinulosa.ALSI',
                           'Dicksonia.antarctica.DIAN', 'Thyrsopteris.elegans.EWXK','Lygodium.flexuosum.LYFL',
                           'Sticherus.truncatus.STTR', 'Osmunda.japonica.OSJA') %>% 
  filter(Azolla.pinnata.AZPN < 100 & Alsophila.spinulosa.ALSI < 100 & 
           Dicksonia.antarctica.DIAN < 100 & Thyrsopteris.elegans.EWXK < 100 &
           Lygodium.flexuosum.LYFL < 100 & Sticherus.truncatus.STTR < 100 & Osmunda.japonica.OSJA < 100)  %>% 
  filter(Osmunda.japonica.OSJA >= 1) %>% 
  filter(Azolla.pinnata.AZPN >=1 | Alsophila.spinulosa.ALSI >=1 |
           Dicksonia.antarctica.DIAN >=1 | Thyrsopteris.elegans.EWXK >=1 | Lygodium.flexuosum.LYFL >=1 |
           Sticherus.truncatus.STTR >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Azolla.pinnata.AZPN, Alsophila.spinulosa.ALSI,
                        Dicksonia.antarctica.DIAN, Thyrsopteris.elegans.EWXK, Lygodium.flexuosum.LYFL, 
                        Sticherus.truncatus.STTR, Osmunda.japonica.OSJA))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.3159669

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("CYATH1_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.3159669,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001373651
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001130859
#### CYATH 2 ####
taxa <- c("Plagiogyria_stenoptera_PLST", "Alsophila_spinulosa_ALSI", "Thyrsopteris_elegans_EWXK", 
          "Lygodium_flexuosum_LYFL", "Sticherus_truncatus_STTR", "Osmunda_japonica_OSJA", "Culcita_macrocarpa_PNZO")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "CYATH2_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("CYATH2_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Plagiogyria.stenoptera.PLST', 'Alsophila.spinulosa.ALSI',
                           'Culcita.macrocarpa.PNZO', 'Thyrsopteris.elegans.EWXK','Lygodium.flexuosum.LYFL',
                           'Sticherus.truncatus.STTR', 'Osmunda.japonica.OSJA') %>% 
  filter(Plagiogyria.stenoptera.PLST < 100 & Alsophila.spinulosa.ALSI < 100 & 
           Culcita.macrocarpa.PNZO < 100 & Thyrsopteris.elegans.EWXK < 100 &
           Lygodium.flexuosum.LYFL < 100 & Sticherus.truncatus.STTR < 100 & Osmunda.japonica.OSJA < 100)  %>% 
  filter(Osmunda.japonica.OSJA >= 1) %>% 
  filter(Plagiogyria.stenoptera.PLST >=1 | Alsophila.spinulosa.ALSI >=1 |
           Culcita.macrocarpa.PNZO >=1 | Thyrsopteris.elegans.EWXK >=1 | Lygodium.flexuosum.LYFL >=1 |
           Sticherus.truncatus.STTR >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Plagiogyria.stenoptera.PLST, Alsophila.spinulosa.ALSI,
                        Culcita.macrocarpa.PNZO, Thyrsopteris.elegans.EWXK, Lygodium.flexuosum.LYFL, 
                        Sticherus.truncatus.STTR, Osmunda.japonica.OSJA))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.302914

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("CYATH2_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.302914,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001433139
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001520378

####CYATH3 ####
taxa <- c("Alsophila_spinulosa_ALSI", "Thyrsopteris_elegans_EWXK", 
          "Lygodium_flexuosum_LYFL", "Dicksonia_antarctica_DIAN", "Cibotium_barometz_CIBA", "Azolla_pinnata_AZPN")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "CYATH3_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("CYATH3_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Alsophila.spinulosa.ALSI',
                           'Dicksonia.antarctica.DIAN', 'Thyrsopteris.elegans.EWXK','Lygodium.flexuosum.LYFL',
                           'Cibotium.barometz.CIBA', 'Azolla.pinnata.AZPN') %>% 
  filter(Alsophila.spinulosa.ALSI < 100 & 
           Dicksonia.antarctica.DIAN < 100 & Thyrsopteris.elegans.EWXK < 100 &
           Lygodium.flexuosum.LYFL < 100 & Cibotium.barometz.CIBA < 100 & Azolla.pinnata.AZPN < 100)  %>% 
  filter(Lygodium.flexuosum.LYFL >= 1) %>% 
  filter(Alsophila.spinulosa.ALSI >=1 |
           Dicksonia.antarctica.DIAN >=1 | Thyrsopteris.elegans.EWXK >=1 | Cibotium.barometz.CIBA >=1 |
           Azolla.pinnata.AZPN >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Alsophila.spinulosa.ALSI,
                        Dicksonia.antarctica.DIAN, Thyrsopteris.elegans.EWXK, Lygodium.flexuosum.LYFL, 
                        Cibotium.barometz.CIBA, Azolla.pinnata.AZPN))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.310327

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("CYATH3_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.310327,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001472288
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001412466

#### LINDS ####
taxa <- c("Lindsaea_microphylla_YIXP", "Odontosira_chinensis_ODCH", "Osmolindsaea_odorata_OSOD", 
          "Lonchitis_hirsuta_VVRN", "Hypolepis_punctata_HYPU", "Alsophila_spinulosa_ALSI", "Azolla_pinnata_AZPN")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "LINDS_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("LINDS_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Lindsaea.microphylla.YIXP', 'Odontosira.chinensis.ODCH',
                           'Osmolindsaea.odorata.OSOD', 'Lonchitis.hirsuta.VVRN','Hypolepis.punctata.HYPU',
                           'Alsophila.spinulosa.ALSI', 'Azolla.pinnata.AZPN') %>% 
  filter(Lindsaea.microphylla.YIXP < 100 & Odontosira.chinensis.ODCH < 100 & 
           Osmolindsaea.odorata.OSOD < 100 & Lonchitis.hirsuta.VVRN < 100 &
           Hypolepis.punctata.HYPU < 100 & Alsophila.spinulosa.ALSI < 100 & Azolla.pinnata.AZPN < 100)  %>% 
  filter(Azolla.pinnata.AZPN >= 1) %>% 
  filter(Lindsaea.microphylla.YIXP >=1 | Odontosira.chinensis.ODCH >=1 |
           Osmolindsaea.odorata.OSOD >=1 | Lonchitis.hirsuta.VVRN >=1 | Hypolepis.punctata.HYPU >=1 |
           Alsophila.spinulosa.ALSI >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Lindsaea.microphylla.YIXP, Odontosira.chinensis.ODCH,
                        Osmolindsaea.odorata.OSOD, Lonchitis.hirsuta.VVRN, Hypolepis.punctata.HYPU, 
                        Alsophila.spinulosa.ALSI, Azolla.pinnata.AZPN))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.215149

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("LINDS_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.215149,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.0021524
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.00200688

#### PTER ####
taxa <- c("Vittaria_appalachiana_NDUV", "Vittaria_lineata_SKYV", 
          "Antrophyum_callifolium_ANCA", "Haplopteris_heterophylla_HAHE", "Adiantum_caudatum_ADCA",
          "Myriopteris_rufa_GSXD", "Pteris_vittata_PTVI")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "PTER_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("PTER_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Vittaria.appalachiania.NDUV', 'Vittaria.lineata.SKYV',
                           'Antrophyum.callifolium.ANCA', 'Haplopteris.heterophylla.HAHE','Adiantum.caudatum.ADCA',
                           'Myriopteris.rufa.GSXD', 'Pteris.vittata.PTVI') %>% 
  filter(Vittaria.appalachiania.NDUV < 100 & Vittaria.lineata.SKYV < 100 &
           Antrophyum.callifolium.ANCA < 100 & Haplopteris.heterophylla.HAHE < 100 &
           Adiantum.caudatum.ADCA < 100 & Myriopteris.rufa.GSXD < 100 & Pteris.vittata.PTVI < 100)  %>% 
  filter(Pteris.vittata.PTVI >= 1) %>% 
  filter(Vittaria.appalachiania.NDUV >=1 |
           Vittaria.lineata.SKYV >=1 | Antrophyum.callifolium.ANCA >=1 | Haplopteris.heterophylla.HAHE >=1 |
           Adiantum.caudatum.ADCA >= 1 | Myriopteris.rufa.GSXD >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Vittaria.appalachiania.NDUV, Vittaria.lineata.SKYV, 
                        Antrophyum.callifolium.ANCA, Haplopteris.heterophylla.HAHE, Adiantum.caudatum.ADCA, 
                        Myriopteris.rufa.GSXD, Pteris.vittata.PTVI))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.02542

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("PTER_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.02542,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.003973515
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.004321253

#### LEPTO ####
taxa <- c("Angiopteris_fokiensis_ANFK", "Equisetum_arvense_EQAR", 
          "Lygodium_japonicum_LYJA", "Osmunda_sp_UOMY_gametophyte", "Psilotum_nudum_PSNU",
          "Salvinia_natans_SANA", "Trichomanes_badium_TRBA")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "LEPTO_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("LEPTO_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Angiopteris.fokiensis.ANFK', 'Equisetum.arvense.EQAR',
                           'Lygodium.japonicum.LYJA', 'Osmunda.sp.UOMY','Psilotum.nudum.PSNU',
                           'Salvinia.natans.SANA', 'Trichomanes.badium.TRBA') %>% 
  filter(Angiopteris.fokiensis.ANFK < 100 & Equisetum.arvense.EQAR < 100 &
           Lygodium.japonicum.LYJA < 100 & Osmunda.sp.UOMY < 100 &
           Psilotum.nudum.PSNU < 100 & Salvinia.natans.SANA < 100 & Trichomanes.badium.TRBA < 100)  %>% 
  filter(Equisetum.arvense.EQAR >= 1) %>% 
  filter(Angiopteris.fokiensis.ANFK >=1 |
           Osmunda.sp.UOMY >=1 | Lygodium.japonicum.LYJA >=1 | Psilotum.nudum.PSNU >=1 |
           Salvinia.natans.SANA >= 1 | Trichomanes.badium.TRBA >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Angiopteris.fokiensis.ANFK, Equisetum.arvense.EQAR, 
                        Lygodium.japonicum.LYJA, Osmunda.sp.UOMY, Psilotum.nudum.PSNU, 
                        Salvinia.natans.SANA, Trichomanes.badium.TRBA))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.232988

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("LEPTO_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.232988,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001036689
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.0007420914

#### POLY ####
taxa <- c("Dryopteris_decipiens_DRDE", "Athyrium_filix_femina_URCP", 
          "Pteridium_aquilinum_subsp._latisculum_PTAQ", "Acrostichum_aureum_ACAR",
          "Lindsaea_linearis_NOKI", "Alsophila_spinulosa_ALSI", "Azolla_pinnata_AZPN", 
          "Lygodium_flexuosum_LYFL")
taxa_to_keep <- data.frame(taxa, row.names = taxa)
name_list <- name.check(time_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(time_tree, names_to_drop)
trimmed_tree <- ladderize(trimmed_tree, right = F)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "POLY_trimmed_tree_toedit.tre")
#Edit tree manually to convert to SIMMAP format 

OrthoFinder <- read.delim("POLY_Orthogroups.GeneCount.tsv")

#Select only taxa of interest 
for_WGDgc <- dplyr::select(OrthoFinder,
                           'Acrostichum.aureum.ACAR', 'Alsophila.spinulosa.ALSI',
                           'Athyrium.filix.femina.URCP', 'Azolla.pinnata.AZPN','Dryopteris.decipiens.DRDE',
                           'Lindsaea.linearis.NOKI', 'Lygodium.flexuosum.LYFL', 
                           'Pteridium.aquilinum.subsp.latisculum.PTAQ') %>% 
  filter(Acrostichum.aureum.ACAR < 100 & Alsophila.spinulosa.ALSI < 100 &
           Athyrium.filix.femina.URCP < 100 & Azolla.pinnata.AZPN < 100 &
           Dryopteris.decipiens.DRDE < 100 & Lindsaea.linearis.NOKI < 100 &
           Lygodium.flexuosum.LYFL < 100 & Pteridium.aquilinum.subsp.latisculum.PTAQ < 100)  %>% 
  filter(Lygodium.flexuosum.LYFL >= 1) %>% 
  filter(Acrostichum.aureum.ACAR >=1 | Alsophila.spinulosa.ALSI >=1 |
           Athyrium.filix.femina.URCP >=1 | Azolla.pinnata.AZPN >=1 | Dryopteris.decipiens.DRDE >=1 |
           Lindsaea.linearis.NOKI >= 1 | Pteridium.aquilinum.subsp.latisculum.PTAQ >= 1)

tmp <- for_WGDgc[rowSums(for_WGDgc[])>0,]

df <- tmp %>% rowwise() %>% 
  mutate(average=mean(c(Lygodium.flexuosum.LYFL, Acrostichum.aureum.ACAR, 
                        Alsophila.spinulosa.ALSI, Athyrium.filix.femina.URCP, Azolla.pinnata.AZPN, 
                        Dryopteris.decipiens.DRDE, Lindsaea.linearis.NOKI, Pteridium.aquilinum.subsp.latisculum.PTAQ))) 

ggplot(data = df, mapping = aes(x = average)) + xlim(0, 3) + geom_density() + theme_classic()

gm_mean(df$average) # 1.190831

#Read in ultrametric species tree in Simmap format 
tree <- phyext::read.simmap("POLY_trimmed_tree_edited.tre", vers = 1.1)
plot(tree)
#Generate subsets and estimate lambda and mu with WGDgc 
birth = vector(mode = "list", length = 10)
death = vector(mode = "list", length = 10)
for (i in 1:10) {
  subset <- sample_n(tmp, 500)
  MLE <- MLEGeneCount(tree, geneCountData = subset, geomMean = 1.190831,
                      conditioning = "oneInBothClades", fixedRetentionRates = T)
  birth[i] <- MLE$birthrate
  death[i] <- MLE$deathrate
}

birth_df <- melt(as.data.frame(birth))
mean(birth_df$value) # 0.001552232 
death_df <- melt(as.data.frame(death))
mean(death_df$value) # 0.001317193

##### Plot MAPS ####

all_maps <- read.csv("MAPS_all.csv")

ggplot(data = all_maps, mapping = aes(x = Node, y = Perc.Dup, color = Type, shape = Type)) + geom_line(size =0.7)+  
  geom_point(size =2) + facet_wrap(~Event, scales = "free_x") + theme_classic() +
  ylab("Percent Duplicated Subtrees") + scale_color_manual(values = c("black","blue", "red")) +
  ylim(0,60)
