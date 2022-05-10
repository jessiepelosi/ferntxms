###########################
## gene_retention.R 
## Jessie Pelosi
## Last modified May 10, 2022
###########################

library(MASS)
library(dplyr)
library(gplots)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(pals)
library(ggExtra)
library(shiny)
library(tibble)
library(reshape2)

#### Plot mean number of duplications by order ####

wgds <- read.csv("WGDperTaxa.csv")

wgds$Order <- factor(wgds$Order, 
                     levels = rev(c("Polypodiales", "Cyatheales",
                                            "Salviniales", "Schizaeales", 
                                            "Hymenophyllales", "Gleicheniales",
                                            "Osmundales","Marattiales", 
                                            "Psilotales", "Ophioglossales", 
                                            "Equisetales")))

ggplot(data = wgds, mapping = aes(y= Order, x = No.Inferred.WGDs, fill= Order)) +
  geom_violin() + xlim(1,3) + theme_classic() +
  stat_summary(fun.y=mean, geom="point", size=2, fill = "black")+
  xlab("Mean Number of Inferred WGDs") + theme(legend.position = "none")

#### Look at copy number distribution? ####

OGs <- read.delim("../OrthoFinder_8.2.21_new/comp_genomics_stats/Orthogroups.GeneCount.tsv",
                  na.strings = 0, stringsAsFactors = F)

OGs_test <- select(OGs,-Amborella_trichopoda.AMTR1.0.pep.all, 
                   -Ginkgo_biloba.HiC, -Arabidopsis_thaliana.TAIR10.pep.all, 
                   -Physcomitrium_patens.Phypa_V3.pep.all, Selaginella_moellendorffii.v1.0.pep.all, 
                   -Total) %>% 
  mutate(number = rowSums(is.na(OGs))) %>% 
  filter(number < 60)

OGs_test_2 <- OGs_test %>% 
  select(-number)
OGs_test_2 <- OGs_test_2[rowSums(OGs_test_2 > 4, na.rm = T) < 2,]
OGs_test_2$k3 <- k3$cluster

OGs_families_to_use <- OGs_test_2 %>% 
  select(-Orthogroup)

small <- OGs_families_to_use[rowSums(OGs_families_to_use > 4, na.rm = T) < 1,]

#big <- rowSums(OGs_families_to_use >= 50, na.rm = TRUE); big <- as.data.frame(big)
#small <- ifelse(big$big > 0, FALSE, TRUE); OGs_families_to_use[,small]
SC <- rowSums(small == 1, na.rm = TRUE); SC <- as.data.frame(SC)
tot <- rowSums(is.na(small)); tot <- as.data.frame(tot)
prop <- SC$SC/(239- tot$tot); prop <- as.data.frame(prop)

ggplot(data = prop, mapping = aes(x = prop, ..scaled..), color = "black") + 
  geom_density(fill = "gray") + theme_classic() +
  xlab("Proportion of Taxa Single Copy") + ylab("Number of Gene Families")

k2 <- kmeans(prop, centers = 2)
k3 <- kmeans(prop, centers = 3)

prop$k3 <- k3$cluster
ggplot() + 
  geom_density(data = prop, mapping = aes(x = prop, y = ..scaled..), fill= "gray", alpha = 0.75) +
  #geom_density(data = prop, mapping = aes(x = prop, fill= as.factor(k3), y=..count..), alpha = 0.5) +
  xlab("Proportion of Taxa Single Copy") + 
  ylab("Gene Family Density") + theme_classic() 

OGs_melt <- melt(OGs_test_2, id.vars = c("k3", "Orthogroup")) %>% 
  mutate(Color = ifelse(value == 2, "blue", ifelse(value < 2, "yellow", ifelse(value >2, "green", "red"))))

#key_table <- read.csv("orders_samples.csv")
#combined <- merge(OGs_melt, key_table, by="variable")
#combined_2 <- combined %>% 
#  group_by(Order, Orthogroup) %>% 
#  summarize(avg = round(mean(value, na.rm = T))) %>% 
#  mutate(Color = ifelse(avg == 2, "blue", ifelse(avg < 2, "yellow", ifelse(avg >2, "green", "red"))))

clrs <- c("green" = "darkgreen", "yellow" = "yellow2", "red" = "red3", "blue" = "dodgerblue2")

OGs_melt$Orthogroup = reorder(OGs_melt$Orthogroup, OGs_melt$k3)

#OGs_melt$Orthogroup <- factor(OGs_melt$Orthogroup, levels=c("Y", "X", "Z"))

ggplot(data = OGs_melt, mapping = aes(x = Orthogroup, y = variable, fill = Color)) +
  geom_bin_2d() + scale_fill_manual(values = clrs)

#### Read in outputs from GOGetter pipeline ####
setwd("C:/Users/Owner/Dropbox/(Insert cool lab name here)/Jessie/Fern_transcriptomes/fern_txms/all_tables/Run_ATH2021_12.4.21/RawCount/")
temp = list.files(pattern="*.cds")
files = lapply(temp, read.csv, header = T)
combined.df <- do.call(cbind, files)
write.csv(file = "ATH2021_all_count_files.csv", x = combined.df, quote = F)

setwd("C:/Users/Owner/Dropbox/(Insert cool lab name here)/Jessie/Fern_transcriptomes/fern_txms/all_tables/Run_ATH2021_12.4.21/Frequency/")
temp = list.files(pattern="*.cds")
files = lapply(temp, read.csv, header = T)
combined.df <- do.call(cbind, files)
write.csv(file = "ATH2021_all_freq_files.csv", x = combined.df, quote = F)

#### Visualizing gene retention w/ heatmap from GOGetter ####
## THIS NEEDS TO BE DONE ON FREQUENCY DATA, NOT RAW COUNTS!!!!!!
## Written by Emily Sessa, rewritten by Jessie Pelosi 

heatmap_data <- read.csv("ATH2021_all_freq_files_EDIT.csv", header = T, row.names = 1)
ordered <- heatmap_data[order(heatmap_data$average, decreasing = T),]
heat_matrix <- as.matrix(ordered)
 
# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_75 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 75, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_75",
  ns = asNamespace("pheatmap")
)

quantile_breaks <- function(xs, n = 96) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(heat_matrix, n = 97)

pheatmap(heat_matrix, cellwidth = 2, cellheight = 2, fontsize = 1, cluster_rows = F, cluster_cols = F, 
        border_color = "black", breaks = mat_breaks, col = tol.rainbow(n=97))

#### PCA on gene retention patterns ####

GO_PCA_data <-  read.csv("ATH2021_all_freq_files_TaxaRows.csv")

PCA_test <- prcomp(GO_PCA_data[,3:97])

autoplot(PCA_test, data = GO_PCA_data, colour = "Type", frame = F, )

PCA <- as.data.frame(PCA_test$x)
groups <- GO_PCA_data$Type
PCA$Type <- groups

scatter <- ggplot(data = PCA, mapping = aes(x=PC1, y=PC2, color = Type)) + geom_point() +
  theme_classic() + xlab("PC1 (52.14%)") + ylab("PC2 (12.94%)") + theme(legend.position =  "none")

ggExtra::ggMarginal(scatter, type = "histogram", groupFill = T)

ggsave("ATH2021_GOPCA.png", dpi = 300, height = 7, width = 7)

#### Compute Linear Discriminant Analysis ####

ind <- sample(2, nrow(GO_PCA_data),
              replace = TRUE,
              prob = c(0.5, 0.5))
training <- GO_PCA_data[ind==1,]
testing <- GO_PCA_data[ind==2,]

lda_training <- training[,2:97]
lda_testing <- testing[,2:97]

linear <- lda(Type~. , lda_training)
linear$prior
linear
p <- predict(linear, lda_testing)
p
ldahist(data = p$x[,1], g = GO_PCA_data$Type)

#### Chi2 Tests ####

counts <- read.csv("../RawCount/ATH2021_Counts_TaxaRows.csv", stringsAsFactors = F)

# Remove txms w/o Ks peaks 
counts_filt <- counts %>% 
  dplyr::filter(Sequence != "Actinostachys_digitata_ACDI" & Sequence != "Aglaomorpha_fortunei_AGFO" &
                  Sequence != "Aleuritopteris_leptolepis_ALLE" & Sequence != "Anemia_phyllitidis_ANPH" &
                  Sequence != "Anemia_tomentosa_CQPW" & Sequence != "Arthropteris_palisotii_ARPA" & 
                  Sequence != "Azolla_caroliniana_CVEG" & Sequence != "Azolla_pinnata_AZPI" & 
                  Sequence != "Azolla_pinnata_AZPN" & Sequence != "Bolbitis_heteroclita_BOHE" & 
                  Sequence != "Davallia_fejeensis_QQWW" & Sequence != "Diplazium_viridescens_DIVI" &
                  Sequence != "Diplazium_esculentum_DIES" & Sequence != "Elaphoglossum_mccluri_ELMC" &
                  Sequence != "Elaphoglossum_yoshinagae_ELYO" & Sequence != "Grammitis_dorsipila_GRDO" &
                  Sequence != "Gymnocarpium_oyamense_GYOY" & Sequence != "Lepisorus_albertii_LEAL" &
                  Sequence != "Lomariopsis_boninensis_LOBO" & Sequence != "Lomariopsis_spectabilis_LOSP" & 
                  Sequence != "Lygodium_flexuosum_LYFL" & Sequence != "Marsilea_quadrifolia_MAQF" & 
                  Sequence != "Marsilea_quadrifolia_MAQU" & Sequence != "Microlepia_hookeriana_MIHO" & 
                  Sequence != "Microlepia_platyphylla_MIPL" & Sequence != "Microlepia_speluncae_MISP" & 
                  Sequence != "Monachosorum_henryi_MOHE" & Sequence != "Monachosorum_maximowiczii_MOMA" & 
                  Sequence != "Oleandra_sp._XQ.2018_OLSP" & Sequence != "Ophioglossum_vulgatum_QHVS" & 
                  Sequence != "Oreogrammitis_congener_ORCO" & Sequence != "Phegopteris_ducursive.pinnata_PHDP" &
                  Sequence != "Phlebodium_pseudoaureum_ZQYU" & Sequence != "Pilularia_globulifera_KIIX" & 
                  Sequence != "Pityrogramma_trifoliata_UJTT" & Sequence != "Plagiogyria_japonica_UWOD" &
                  Sequence != "Polypodium_hesperium_GYFU" & Sequence != "Polypodium_hesperium_IXLH" & 
                  Sequence != "Polypodium_virginianum_POVI" & Sequence != "Prosaptia_obliquata_PROB" &
                  Sequence != "Pteris_fauriei_PTFA" & Sequence != "Pteris_vittata_POPJ" & Sequence != "Pteris_vittata_PTVI" &
                  Sequence != "Salvinia_natans_SANA" & Sequence != "Salvinia_natans_SANT" & Sequence != "Stenogramma_wilfordii_STNO" &
                  Sequence != "Vandenboschia_striata_VAST" & Sequence != "Woodsia_ilvensis_WOIL" &
                  Sequence != "Didymochlaena_runcatula_RFRB" & Sequence != "Lomagramma_sumatrana_LOSU" &
                  Sequence != "Ophioglossum_vulgatum_OPVU")

lst <- group_by(counts_filt, Taxon) %>% 
  group_split() 

# Run Chi2 tests on each group comparing GO compositions of peaks and full transcriptomes 
for (i in 1:length(lst)) {
  tmp_df1 <- lst[[i]]
  tmp_df2 <- dplyr::select(tmp_df1, !(Taxon:Type)) %>% 
    tibble::column_to_rownames(var = "Sequence")
  c2t <- chisq.test(as.matrix(tmp_df2))
  df <- as.data.frame(c2t$residuals) %>% 
    rownames_to_column() %>% rbind(df)
  print(c2t$p.value)
  print(c2t$statistic)
}

write.csv(df, file = "single_peak_resids.csv", quote = F)

residuals <- read.csv("all_tables/Run_ATH2021_12.4.21/Frequency/paralog_resids.csv", header = T)

residuals_melt <- melt(residuals)

residuals_long <- melt(residuals) %>% 
  mutate(average = mean(value)) %>% 
  mutate(Color = ifelse(value > 2, "red", ifelse(value < -2, "blue", "gray"))) %>% 
  filter(rowname != "Mean")

cols <- c("red" = "red", "blue"= "dodgerblue2", "gray"= "lightgray")

level_order <- as.vector(avg$variable)

ggplot(data = residuals_long, mapping = aes(x = rowname, y = factor(variable, level = level_order), fill = Color)) + geom_bin_2d() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Plot by event 
residuals_event <- mutate(residuals_melt, average = mean(value)) %>% 
  group_by(event, variable) %>% 
  summarize(average=mean(value)) %>% 
  mutate(Color = ifelse(average > 2, "red", ifelse(average < -2, "blue", "gray"))) %>% 
  filter(event == "OPHIO.1" | event == "DIPT.2" | event == "PSIL.2" | event =="CYATH.2" |
         event == "CYATH.3" | event ==  "PTER.3" | event == "LINDS" | event == "EQUI" | 
         event == "GLEI" | event == "MARA" | event == "DIPT.1" | event == "OPHIO.2" |
        event == "PSIL.1" | event == "HYMN" | event ==  "CYATH.1" |
         event == "POLY" | event == "LEPTO" | event == "SALV")

residuals_event$event <- factor(residuals_event$event,levels = c("OPHIO.1", "DIPT.2", "PSIL.2", "CYATH.2", "CYATH.3", "PTER.3", "LINDS", "EQUI", "GLEI", "MARA", "DIPT.1", "OPHIO.2", "PSIL.1", "HYMN", "CYATH.1", "POLY", "LEPTO", "SALV"))


# Sort by average value 

residuals$event

avg <- residuals %>% 
  filter(event == "Mean") %>% 
  melt() %>% 
  arrange(value) %>% 
  select(variable)

level_order <- as.vector(avg$variable)

ggplot(data = residuals_event, mapping = aes(x = event, y = factor(variable, level = level_order), fill = Color)) +
  geom_bin_2d() + scale_fill_manual(values = cols) + 
  ylab(label = "GO Slim Category") + xlab("WGD Event") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=1)) +
  theme(axis.text.y = element_text(size = 7))


# For mutliple peaks- do individually :( 

##### Adiantum raddianum BMJR ####
F1 <- counts %>% 
  filter(Taxon == "Adiantum_raddianum_BMJR") %>% 
  filter(Sequence != "Adiantum_raddianum_BMJR_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Adiantum_raddianum_BMJR_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Adiantum_raddianum_BMJR") %>% 
  filter(Sequence != "Adiantum_raddianum_BMJR_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Adiantum_raddianum_BMJR_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Adiantum_raddianum_BMJR") %>% 
  filter(Sequence != "Adiantum_raddianum_BMJR") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test

##### Alsophila podophylla ASPD #####
F1 <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPD") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPD_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_podophylla_ASPD_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPD") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPD_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_podophylla_ASPD_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPD") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPD") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Alsophila podophylla ASPO ####
F1 <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPO") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPO_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_podophylla_ASPO_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPO") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPO_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_podophylla_ASPO_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Alsophila_podophylla_ASPO") %>% 
  filter(Sequence != "Alsophila_podophylla_ASPO") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Alsophila_sp._XQ.2018_ALSP ####
F1 <- counts %>% 
  filter(Taxon == "Alsophila_sp._XQ.2018_ALSP") %>% 
  filter(Sequence != "Alsophila_sp_ALSP_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_sp_ALSP_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Alsophila_sp._XQ.2018_ALSP") %>% 
  filter(Sequence != "Alsophila_sp_ALSP_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_sp_ALSP_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Alsophila_sp._XQ.2018_ALSP") %>% 
  filter(Sequence != "Alsophila_sp._XQ.2018_ALSP") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Alsophila spinulosa ALSI ####
F1 <- counts %>% 
  filter(Taxon == "Alsophila_spinulosa_ALSI") %>% 
  filter(Sequence != "Alsophila_spinulosa_ALSI_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_spinulosa_ALSI_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Alsophila_spinulosa_ALSI") %>% 
  filter(Sequence != "Alsophila_spinulosa_ALSI_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Alsophila_spinulosa_ALSI_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Alsophila_spinulosa_ALSI") %>% 
  filter(Sequence != "Alsophila_spinulosa_ALSI") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Cibotium barometz CIBA #####
F1 <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBA") %>% 
  filter(Sequence != "Cibotium_barometz_CIBA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cibotium_barometz_CIBA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBA") %>% 
  filter(Sequence != "Cibotium_barometz_CIBA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cibotium_barometz_CIBA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBA") %>% 
  filter(Sequence != "Cibotium_barometz_CIBA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Cibotium barometz CIBZ #####
F1 <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBZ") %>% 
  filter(Sequence != "Cibotium_barometz_CIBZ_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cibotium_barometz_CIBZ_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBZ") %>% 
  filter(Sequence != "Cibotium_barometz_CIBZ_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cibotium_barometz_CIBZ_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Cibotium_barometz_CIBZ") %>% 
  filter(Sequence != "Cibotium_barometz_CIBZ") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Cyathea spinulosa GANB ####
F1 <- counts %>% 
  filter(Taxon == "Cyathea_spinulosa_GANB") %>% 
  filter(Sequence != "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Cyathea_spinulosa_GANB") %>% 
  filter(Sequence != "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Cyathea_spinulosa_GANB") %>% 
  filter(Sequence != "Cyathea_spinulosa_GANB") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Dicksonia antarctica DIAN #####
F1 <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAN") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAN_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dicksonia_antarctica_DIAN_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAN") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAN_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dicksonia_antarctica_DIAN_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAN") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAN") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Dicksonia antarctica DIAT #####
F1 <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAT") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAT_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dicksonia_antarctica_DIAT_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAT") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAT_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dicksonia_antarctica_DIAT_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Dicksonia_antarctica_DIAT") %>% 
  filter(Sequence != "Dicksonia_antarctica_DIAT") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Diploptergium laevissimum DILA ####
F1 <- counts %>% 
  filter(Taxon == "Diplopterygium_laevissimum_DILA") %>% 
  filter(Sequence != "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Diplopterygium_laevissimum_DILA") %>% 
  filter(Sequence != "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Diplopterygium_laevissimum_DILA") %>% 
  filter(Sequence != "Diplopterygium_laevissimum_DILA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Dipteris lobbiana DILO ####
F1 <- counts %>% 
  filter(Taxon == "Dipteris_lobbiana_DILO") %>% 
  filter(Sequence != "Dipteris_lobiana_DILO_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dipteris_lobiana_DILO_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Dipteris_lobbiana_DILO") %>% 
  filter(Sequence != "Dipteris_lobiana_DILO_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Dipteris_lobiana_DILO_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Dipteris_lobbiana_DILO") %>% 
  filter(Sequence != "Dipteris_lobbiana_DILO") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Haplopteris amboinensis HAAM #####
F1 <- counts %>% 
  filter(Taxon == "Haplopteris_amboinensis_HAAM") %>% 
  filter(Sequence != "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Haplopteris_amboinensis_HAAM") %>% 
  filter(Sequence != "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Haplopteris_amboinensis_HAAM") %>% 
  filter(Sequence != "Haplopteris_amboinensis_HAAM") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Haplopteris heterophylla HAHE ####
F1 <- counts %>% 
  filter(Taxon == "Haplopteris_heterophylla_HAHE") %>% 
  filter(Sequence != "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Haplopteris_heterophylla_HAHE") %>% 
  filter(Sequence != "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Haplopteris_heterophylla_HAHE") %>% 
  filter(Sequence != "Haplopteris_heterophylla_HAHE") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Lindsea heterophylla LIHE ####
F1 <- counts %>% 
  filter(Taxon == "Lindsaea_heterophylla_LIHE") %>% 
  filter(Sequence != "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Lindsaea_heterophylla_LIHE") %>% 
  filter(Sequence != "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Lindsaea_heterophylla_LIHE") %>% 
  filter(Sequence != "Lindsaea_heterophylla_LIHE") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Lindsea linearis NOKI #####
F1 <- counts %>% 
  filter(Taxon == "Lindsaea_linearis_NOKI") %>% 
  filter(Sequence != "Lindsaea_linearis_NOKI_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_linearis_NOKI_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Lindsaea_linearis_NOKI") %>% 
  filter(Sequence != "Lindsaea_linearis_NOKI_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_linearis_NOKI_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Lindsaea_linearis_NOKI") %>% 
  filter(Sequence != "Lindsaea_linearis_NOKI") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Lindsea microphylla YIXP #####
F1 <- counts %>% 
  filter(Taxon == "Lindsaea_microphylla_YIXP") %>% 
  filter(Sequence != "Lindsaea_microphylla_YIXP_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_microphylla_YIXP_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Lindsaea_microphylla_YIXP") %>% 
  filter(Sequence != "Lindsaea_microphylla_YIXP_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lindsaea_microphylla_YIXP_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Lindsaea_microphylla_YIXP") %>% 
  filter(Sequence != "Lindsaea_microphylla_YIXP") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Lonchitis hirsuta VVRN #####
F1 <- counts %>% 
  filter(Taxon == "Lonchitis_hirsuta_VVRN") %>% 
  filter(Sequence != "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Lonchitis_hirsuta_VVRN") %>% 
  filter(Sequence != "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Lonchitis_hirsuta_VVRN") %>% 
  filter(Sequence != "Lonchitis_hirsuta_VVRN") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Plagiogyria_japonica_PLJA ####
F1 <- counts %>% 
  filter(Taxon == "Plagiogyria_japonica_PLJA") %>% 
  filter(Sequence != "Plagiogyria_japonica_PLJA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Plagiogyria_japonica_PLJA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Plagiogyria_japonica_PLJA") %>% 
  filter(Sequence != "Plagiogyria_japonica_PLJA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Plagiogyria_japonica_PLJA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Plagiogyria_japonica_PLJA") %>% 
  filter(Sequence != "Plagiogyria_japonica_PLJA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Plagiogyria stenoptera PLST #####
F1 <- counts %>% 
  filter(Taxon == "Plagiogyria_stenoptera_PLST") %>% 
  filter(Sequence != "Plagiogyria_stenoptera_PLST_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Plagiogyria_stenoptera_PLST_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Plagiogyria_stenoptera_PLST") %>% 
  filter(Sequence != "Plagiogyria_stenoptera_PLST_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Plagiogyria_stenoptera_PLST_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Plagiogyria_stenoptera_PLST") %>% 
  filter(Sequence != "Plagiogyria_stenoptera_PLST") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Psilotum nudum PSND ####
F1 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSND") %>% 
  filter(Sequence != "Psilotum_nudum_PSND_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_PSND_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSND") %>% 
  filter(Sequence != "Psilotum_nudum_PSND_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_PSND_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSND") %>% 
  filter(Sequence != "Psilotum_nudum_PSND") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Psilotum nudum PSNU #####
F1 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSNU") %>% 
  filter(Sequence != "Psilotum_nudum_PSNU_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_PSNU_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSNU") %>% 
  filter(Sequence != "Psilotum_nudum_PSNU_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_PSNU_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Psilotum_nudum_PSNU") %>% 
  filter(Sequence != "Psilotum_nudum_PSNU") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Psilotum nudum QVMR #####
F1 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_QVMR") %>% 
  filter(Sequence != "Psilotum_nudum_QVMR_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_QVMR_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Psilotum_nudum_QVMR") %>% 
  filter(Sequence != "Psilotum_nudum_QVMR_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Psilotum_nudum_QVMR_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Psilotum_nudum_QVMR") %>% 
  filter(Sequence != "Psilotum_nudum_QVMR") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Sceptridium dissectum EEAQ #####
F1 <- counts %>% 
  filter(Taxon == "Sceptridium_dissectum_EEAQ") %>% 
  filter(Sequence != "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Sceptridium_dissectum_EEAQ") %>% 
  filter(Sequence != "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Sceptridium_dissectum_EEAQ") %>% 
  filter(Sequence != "Sceptridium_dissectum_EEAQ") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Sceptiridum japonicum SCJA ####
F1 <- counts %>% 
  filter(Taxon == "Sceptridium_japonicum_SCJA") %>% 
  filter(Sequence != "Sceptridium_japonicum_SCJA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Sceptridium_japonicum_SCJA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Sceptridium_japonicum_SCJA") %>% 
  filter(Sequence != "Sceptridium_japonicum_SCJA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Sceptridium_japonicum_SCJA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Sceptridium_japonicum_SCJA") %>% 
  filter(Sequence != "Sceptridium_japonicum_SCJA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Stenochlaena palustris STPA #####
F1 <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPA") %>% 
  filter(Sequence != "Stecnochlaena_palustris_STPA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Stecnochlaena_palustris_STPA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPA") %>% 
  filter(Sequence != "Stecnochlaena_palustris_STPA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Stecnochlaena_palustris_STPA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPA") %>% 
  filter(Sequence != "Stenochlaena_palustris_STPA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Stenochleana palustris STPL #####
F1 <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPL") %>% 
  filter(Sequence != "Stecnochlaena_palustris_STPL_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Stecnochlaena_palustris_STPL_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPL") %>% 
  filter(Sequence != "Stecnochlaena_palustris_STPL_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Stecnochlaena_palustris_STPL_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Stenochlaena_palustris_STPL") %>% 
  filter(Sequence != "Stenochlaena_palustris_STPL") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Taenitis blechnoides TABL ####
F1 <- counts %>% 
  filter(Taxon == "Taenitis_blechnoides_TABL") %>% 
  filter(Sequence != "Taenitis_blechnoides_TABL_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Taenitis_blechnoides_TABL_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Taenitis_blechnoides_TABL") %>% 
  filter(Sequence != "Taenitis_blechnoides_TABL_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Taenitis_blechnoides_TABL_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Taenitis_blechnoides_TABL") %>% 
  filter(Sequence != "Taenitis_blechnoides_TABL") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Tmesipteris tannensis TMTA ####
F1 <- counts %>% 
  filter(Taxon == "Tmesipteris_tannensis_TMTA") %>% 
  filter(Sequence != "Tmesipteris_tannensis_TMTA_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Tmesipteris_tannensis_TMTA_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Tmesipteris_tannensis_TMTA") %>% 
  filter(Sequence != "Tmesipteris_tannensis_TMTA_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Tmesipteris_tannensis_TMTA_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Tmesipteris_tannensis_TMTA") %>% 
  filter(Sequence != "Tmesipteris_tannensis_TMTA") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

##### Vittaria appalachiana NDUV #####
F1 <- counts %>% 
  filter(Taxon == "Vittaria_appalachiana_NDUV") %>% 
  filter(Sequence != "Vittaria_appalachiana_NDUV_WGD_paralogs_peak2") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F1))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Vittaria_appalachiana_NDUV_WGD_paralogs_peak1") %>% 
  rbind(df)

F2 <- counts %>% 
  filter(Taxon == "Vittaria_appalachiana_NDUV") %>% 
  filter(Sequence != "Vittaria_appalachiana_NDUV_WGD_paralogs_peak1") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type))  
test <- chisq.test(as.matrix(F2))
test
test$p.value
df <- as.data.frame(test$residuals) %>% 
  rownames_to_column() %>% filter(rowname == "Vittaria_appalachiana_NDUV_WGD_paralogs_peak2") %>% 
  rbind(df)

P <- counts %>% 
  filter(Taxon == "Vittaria_appalachiana_NDUV") %>% 
  filter(Sequence != "Vittaria_appalachiana_NDUV") %>% 
  column_to_rownames(var = "Sequence") %>% 
  dplyr::select(!(Taxon:Type)) %>% 
  dplyr::select_if(colSums(.) != 0)
test <- chisq.test(as.matrix(P))
test
test$p.value

write.csv(df, file = "mult_peaks_resids.csv", quote = F)
