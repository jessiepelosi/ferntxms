## Ks_plots.R
## Written by Jessie Pelosi for:
## Pelosi, Kim, Barbazuk, and Sessa: Phylotranscriptomics illuminates the placement of whole genome duplications and gene retention in ferns 
## Last modified Septebmer 10, 2021
##

library(ggplot2)
library(dplyr)
library(mixtools)

#### Acrostichum aureum ACAR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.aureum_ACAR <- read.delim("ks_distributions/Acrostichum_aureum_ACAR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.aureum_ACAR_filt <- A.aureum_ACAR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.aureum_ACAR_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.aureum_ACAR_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Acrostichum aureum (ACAR)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.aureum_ACAR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.aureum_ACAR_filt.0 <- A.aureum_ACAR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.aureum_ACAR_normalmixEM <- normalmixEM(A.aureum_ACAR_filt.0$Ks, k=2)
summary(A.aureum_ACAR_normalmixEM)
# summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.124204 0.875796
#mu     0.167943 2.391196
#sigma  0.137978 0.948064
#loglik at estimate:  -6787.756 

# Generate list of paralogs 

A.aureum_ACAR_WGD_paralogs <- A.aureum_ACAR %>% 
  filter(Ks > 1.443132) %>% 
  filter(Ks < 3.33926) %>% 
  select(X)

write.table(A.aureum_ACAR_WGD_paralogs, quote = F, file = "Acrostichum_aureum_ACAR_WGD_paralogs.tsv")

#### Acrostichum aureum ACAU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.aureum_ACAU <- read.delim("ks_distributions/Acrostichum_aureum_ACAU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.aureum_ACAU_filt <- A.aureum_ACAU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.aureum_ACAU_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.aureum_ACAU_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Acrostichum aureum (ACAU)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.aureum_ACAU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.aureum_ACAU_filt.0 <- A.aureum_ACAU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.aureum_ACAU_normalmixEM <- normalmixEM(A.aureum_ACAU_filt.0$Ks, k=2)
summary(A.aureum_ACAU_normalmixEM)
#summary of normalmixEM object:
#        comp 1    comp 2
#lambda 0.896571 0.1034289
#mu     2.383707 0.0962037
#sigma  0.955549 0.0835573
#loglik at estimate:  -5861.032 

# Generate list of paralogs 

A.aureum_ACAU_WGD_paralogs <- A.aureum_ACAU %>% 
  filter(Ks > 1.428158) %>% 
  filter(Ks < 3.339256) %>% 
  select(X)

write.table(A.aureum_ACAU_WGD_paralogs, quote = F, file = "Acrostichum_aureum_ACAU_WGD_paralogs.tsv")

#### Actinostachys digitata ACDI ####
# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.digitata_ACDI <- read.delim("ks_distributions/Actinostachys_digitata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.digitata_ACDI_filt <- A.digitata_ACDI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.digitata_ACDI_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.digitata_ACDI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Actinostachys digitata (ACDI)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.digitata_ACDI.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Acystopteris japonica ACJA #####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.japonica_ACJA <- read.delim("ks_distributions/Acystopteris_japonica.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.japonica_ACJA_filt <- A.japonica_ACJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.japonica_ACJA_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.japonica_ACJA_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Acystopteris japonica (ACJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.japonica_ACJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.japonica_ACJA_filt.0 <- A.japonica_ACJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.japonica_ACJA_normalmixEM <- normalmixEM(A.japonica_ACJA_filt.0$Ks, k=2)
summary(A.japonica_ACJA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.106720 0.893280
#mu     0.133371 2.085389
#sigma  0.109362 0.964866
#loglik at estimate:  -7860.668

# Generate list of paralogs 

A.japonica_ACJA_WGD_paralogs <- A.japonica_ACJA %>% 
  filter(Ks > 1.120523) %>% 
  filter(Ks < 3.050255) %>% 
  select(X)

write.table(A.japonica_ACJA_WGD_paralogs, quote = F, file = "Acystopteris_japonica_ACJA_WGD_paralogs.tsv")

#### Acystopteris tenuisecta ACTE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.tenuisecta_ACTE <- read.delim("ks_distributions/Acystopteris_tenuisecta.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.tenuisecta_ACTE_filt <- A.tenuisecta_ACTE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.tenuisecta_ACTE_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.tenuisecta_ACTE_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Acystopteris tenuisecta (ACTE)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.tenuisecta_ACTE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.tenuisecta_ACTE_filt.0 <- A.tenuisecta_ACTE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.tenuisecta_ACTE_normalmixEM <- normalmixEM(A.tenuisecta_ACTE_filt.0$Ks, k=2)
summary(A.tenuisecta_ACTE_normalmixEM)
#summary of normalmixEM object:
#       comp 1    comp 2
#lambda 0.922728 0.0772717
#mu     2.150663 0.1199264
#sigma  0.958628 0.0988928
#loglik at estimate:  -7152.907  

# Generate list of paralogs 

A.tenuisecta_ACTE_WGD_paralogs <- A.tenuisecta_ACTE %>% 
  filter(Ks > 1.192035) %>% 
  filter(Ks < 3.109291) %>% 
  select(X)

write.table(A.tenuisecta_ACTE_WGD_paralogs, quote = F, file = "Acystopteris_tenuisecta_ACTE_WGD_paralogs.tsv")

#### Adiantum aleucticum WCLG ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.aleuticum_WCLG <- read.delim("ks_distributions/Adiantum_aleuticum_WCLG.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.aleuticum_WCLG_filt <- A.aleuticum_WCLG %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.aleuticum_WCLG_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.aleuticum_WCLG_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Adiantum aleuticum (WCLG)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.aleuticum_WCLG.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.aleuticum_WCLG_filt.0 <- A.aleuticum_WCLG_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.aleuticum_WCLG_normalmixEM <- normalmixEM(A.aleuticum_WCLG_filt.0$Ks, k=2)
summary(A.aleuticum_WCLG_normalmixEM)

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.117995 0.882005
#mu     0.143241 2.158139
#sigma  0.103570 0.975148
#loglik at estimate:  -7059.15 

# Generate list of paralogs 

A.aleuticum_WCLG_WGD_paralogs <- A.aleuticum_WCLG %>% 
  filter(Ks > 1.183) %>% 
  filter(Ks < 3.133) %>% 
  select(X)

write.table(A.aleuticum_WCLG_WGD_paralogs, quote = F, file = "Adiantum_aleuticum_WCLG_WGD_paralogs.tsv")

#### Adiantum capillus-veneris ADVE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.capillus_veneris_ADVE <- read.delim("ks_distributions/Adiantum_capillus-veneris.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.capillus_veneris_ADVE_filt <- A.capillus_veneris_ADVE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.capillus_veneris_ADVE_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.capillus_veneris_ADVE_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Adiantum capillus-veneris (ADVE)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.capillus_veneris_ADVE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.capillus_veneris_ADVE_filt.0 <- A.capillus_veneris_ADVE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.capillus_veneris_ADVE_normalmixEM <- normalmixEM(A.capillus_veneris_ADVE_filt.0$Ks, k=2)
summary(A.capillus_veneris_ADVE_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.106576 0.893424
#mu     0.152587 2.271555
#sigma  0.134767 0.957280
#loglik at estimate:  -7171.028 

# Generate list of paralogs 

A.capillus_veneris_ADVE_WGD_paralogs <- A.capillus_veneris_ADVE %>% 
  filter(Ks > 1.314275) %>% 
  filter(Ks < 3.228835) %>% 
  select(X)

write.table(A.capillus_veneris_ADVE_WGD_paralogs, quote = F, file = "Adiantum_capillus_veneris_ADVE_WGD_paralogs.tsv")

#### Adiantum caudatum ACDA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.caudatum_ACDA <- read.delim("ks_distributions/Adiantum_caudatum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.caudatum_ACDA_filt <- A.caudatum_ACDA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.caudatum_ACDA_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.caudatum_ACDA_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Adiantum caudatum (ACDA)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.caudatum_ACDA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.caudatum_ACDA_filt.0 <- A.caudatum_ACDA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.caudatum_ACDA_normalmixEM <- normalmixEM(A.caudatum_ACDA_filt.0$Ks, k=2)
summary(A.caudatum_ACDA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1105069 0.889493
#mu     0.1001760 2.262920
#sigma  0.0877442 0.969913
#loglik at estimate:  -5632.742

# Generate list of paralogs 

A.caudatum_ACDA_WGD_paralogs <- A.caudatum_ACDA %>% 
  filter(Ks > 1.293007) %>% 
  filter(Ks < 3.232833) %>% 
  select(X)

write.table(A.caudatum_ACDA_WGD_paralogs, quote = F, file = "Adiantum_caudatum_ACDA_WGD_paralogs.tsv")

#### Adiantum raddianum BMJR #### 

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.raddianum_BMJR <- read.delim("ks_distributions/Adiantum_raddianum_BMJR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.raddianum_BMJR_filt <- A.raddianum_BMJR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.raddianum_BMJR_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.raddianum_BMJR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Adiantum raddianum (BMJR)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.raddianum_BMJR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.raddianum_BMJR_filt.0 <- A.raddianum_BMJR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.raddianum_BMJR_normalmixEM <- normalmixEM(A.raddianum_BMJR_filt$Ks, k=2)
# Fits to 2 WGD peaks -- does not fit to recent duplicates 
summary(A.raddianum_BMJR_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.449256 0.550745
#mu     0.474869 2.451718
#sigma  0.297180 0.856945
#loglik at estimate:  -8634.009 

# Generate list of paralogs 

A.raddianum_BMJR_WGD_paralogs_peak1 <- A.raddianum_BMJR %>% 
  filter(Ks > 0.177689) %>% 
  filter(Ks < 0.772049) %>% 
  select(X)

write.table(A.raddianum_BMJR_WGD_paralogs_peak1, quote = F, file = "Adiantum_raddianum_BMJR_WGD_paralogs_peak1.tsv")

A.raddianum_BMJR_WGD_paralogs_peak2 <- A.raddianum_BMJR %>% 
  filter(Ks > 1.594773) %>% 
  filter(Ks < 3.308663) %>% 
  select(X)

write.table(A.raddianum_BMJR_WGD_paralogs_peak2, quote = F, file = "Adiantum_raddianum_BMJR_WGD_paralogs_peak2.tsv")

#### Aglaomorpha bonii AGBO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.bonii_AGBO <- read.delim("ks_distributions/Aglaomorpha_bonii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.bonii_AGBO_filt <- A.bonii_AGBO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.bonii_AGBO_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.bonii_AGBO_filt, mapping = aes(x=Ks, ..scaled..*400), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Aglaomorpha bonii (AGBO)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.bonii_AGBO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.bonii_AGBO_filt.0 <- A.bonii_AGBO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.bonii_AGBO_normalmixEM <- normalmixEM(A.bonii_AGBO_filt.0$Ks, k=2)
# Fits to 2 WGD peaks -- does not fit to recent duplicates 
summary(A.bonii_AGBO_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.153061 0.846939
#mu     0.272712 2.200406
#sigma  0.209759 0.949078
#loglik at estimate:  -5849.769 

# Generate list of paralogs 

A.bonii_AGBOR_WGD_paralogs <- A.bonii_AGBO %>% 
  filter(Ks > 1.251328) %>% 
  filter(Ks < 3.149484) %>% 
  select(X)

write.table(A.bonii_AGBOR_WGD_paralogs, quote = F, file = "Aglaomorpha_bonii_AGBOR_WGD_paralogs_peak1.tsv")

#### Aglaomorpha fortunei AGFO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.fortunei_AGFO <- read.delim("ks_distributions/Aglaomorpha_fortunei.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.fortunei_AGFO_filt <- A.fortunei_AGFO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.fortunei_AGFO_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.fortunei_AGFO_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Aglaomorpha fortunei (AGFO)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.fortunei_AGFO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.fortunei_AGFO_filt.0 <- A.fortunei_AGFO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.fortunei_AGFO_normalmixEM <- normalmixEM(A.fortunei_AGFO_filt.0$Ks, k=2)
summary(A.fortunei_AGFO_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.147280 0.852720
#mu     0.115400 2.148040
#sigma  0.104794 0.998979
#loglik at estimate:  -8763.263

# Generate list of paralogs 

A.fortunei_AGFO_WGD_paralogs <- A.fortunei_AGFO %>% 
  filter(Ks > 1.149061) %>% 
  filter(Ks < 3.147019) %>% 
  select(X)

write.table(A.fortunei_AGFO_WGD_paralogs, quote = F, file = "Aglaomorpha_fortunei_AGFO_WGD_paralogs_peak1.tsv")

#### Aleuritopteris chrysophylla ALCH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.chrysophylla_ALCH <- read.delim("ks_distributions/Aleuritopteris_chrysophylla.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.chrysophylla_ALCH_filt <- A.chrysophylla_ALCH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.chrysophylla_ALCH_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.chrysophylla_ALCH_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Aleuritopteris chrysophylla (ALCH)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.chrysophylla_ALCH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.chrysophylla_ALCH_filt.0 <- A.chrysophylla_ALCH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.chrysophylla_ALCH_normalmixEM <- normalmixEM(A.chrysophylla_ALCH_filt.0$Ks, k=2)
# Fits to 2 WGD peaks -- does not fit to recent duplicates 
summary(A.chrysophylla_ALCH_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.108272 0.891728
#mu     0.130992 2.363125
#sigma  0.115001 0.934159
#loglik at estimate:  -6785.173 

# Generate list of paralogs 

A.chrysophylla_ALCH_WGD_paralogs <- A.chrysophylla_ALCH %>% 
  filter(Ks > 1.428966) %>% 
  filter(Ks < 3.297284) %>% 
  select(X)

write.table(A.chrysophylla_ALCH_WGD_paralogs, quote = F, file = "Aleuritopteris_chrysophylla_ALCH_WGD_paralogs.tsv")

#### Aleuritopteris leptolepis ALLE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.leptolepis_ALLE <- read.delim("ks_distributions/Aleuritopteris_leptolepis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.leptolepis_ALLE_filt <- A.leptolepis_ALLE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.leptolepis_ALLE_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.leptolepis_ALLE_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Aleuritopteris leptolepis (ALLE)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.leptolepis_ALLE.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Alsophila podophylla ASPD ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.podophylla_ASPD <- read.delim("ks_distributions/Alsophila_podophylla_ASPD.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.podophylla_ASPD_filt <- A.podophylla_ASPD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.podophylla_ASPD_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.podophylla_ASPD_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Alsophila podophylla (ASPD)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.podophylla_ASPD.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.podophylla_ASPD_filt.0 <- A.podophylla_ASPD_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.podophylla_ASPD_normalmixEM <- normalmixEM(A.podophylla_ASPD_filt.0$Ks, k=3)
# Fits to 2 WGD peaks -- does not fit to recent duplicates 
summary(A.podophylla_ASPD_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.328309 0.417351 0.254340
#mu     0.316740 1.248150 2.940376
#sigma  0.136046 0.534323 0.593977
#loglik at estimate:  -5246.811 

# Generate list of paralogs 

A.podophylla_ASPD_WGD_paralogs_peak1 <- A.podophylla_ASPD %>% 
  filter(Ks > 0.180694) %>% 
  filter(Ks < 0.452786) %>% 
  select(X)

write.table(A.podophylla_ASPD_WGD_paralogs_peak1, quote = F, file = "Alsophila_podophylla_ASPD_WGD_paralogs_peak1.tsv")

A.podophylla_ASPD_WGD_paralogs_peak2 <- A.podophylla_ASPD %>% 
  filter(Ks > 0.713827) %>% 
  filter(Ks < 1.782473) %>% 
  select(X)

write.table(A.podophylla_ASPD_WGD_paralogs_peak2, quote = F, file = "Alsophila_podophylla_ASPD_WGD_paralogs_peak2.tsv")

#### Alsophila podophylla ASPO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.podophylla_ASPO <- read.delim("ks_distributions/Alsophila_podophylla_ASPO.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.podophylla_ASPO_filt <- A.podophylla_ASPO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.podophylla_ASPO_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.podophylla_ASPO_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Alsophila podophylla (ASPO)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.podophylla_ASPO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.podophylla_ASPO_filt.0 <- A.podophylla_ASPO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.podophylla_ASPO_normalmixEM <- normalmixEM(A.podophylla_ASPO_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3 
summary(A.podophylla_ASPO_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.348356 0.349239 0.302405
#mu     0.295957 1.177717 2.718571
#sigma  0.138891 0.490455 0.693191
#loglik at estimate:  -8661.821 

# Generate list of paralogs 

A.podophylla_ASPO_WGD_paralogs_peak1 <- A.podophylla_ASPO %>% 
  filter(Ks > 0.157066) %>% 
  filter(Ks < 0.434848) %>% 
  select(X)

write.table(A.podophylla_ASPO_WGD_paralogs_peak1, quote = F, file = "Alsophila_podophylla_ASPO_WGD_paralogs_peak1.tsv")

A.podophylla_ASPO_WGD_paralogs_peak2 <- A.podophylla_ASPO %>% 
  filter(Ks > 0.687262) %>% 
  filter(Ks < 1.668172) %>% 
  select(X)

write.table(A.podophylla_ASPO_WGD_paralogs_peak2, quote = F, file = "Alsophila_podophylla_ASPO_WGD_paralogs_peak2.tsv")

#### Alsophila sp. XQ-2018 ALSP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.sp_ALSP <- read.delim("ks_distributions/Alsophila_sp._XQ-2018.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.sp_ALSP_filt <- A.sp_ALSP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.sp_ALSP_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.sp_ALSP_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Alsophila sp. (ALSP)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.sp_ALSP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.sp_ALSP_filt.0 <- A.sp_ALSP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.sp_ALSP_normalmixEM <- normalmixEM(A.sp_ALSP_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3 
summary(A.sp_ALSP_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.351494 0.431694 0.216813
#mu     0.301295 1.351928 3.056959
#sigma  0.148865 0.579019 0.542983
#loglik at estimate:  -9547.274 

# Generate list of paralogs 

A.sp_ALSP_WGD_paralogs_peak1 <- A.sp_ALSP %>% 
  filter(Ks > 0.15243) %>% 
  filter(Ks < 0.45016) %>% 
  select(X)

write.table(A.sp_ALSP_WGD_paralogs_peak1, quote = F, file = "Alsophila_sp_ALSP_WGD_paralogs_peak1.tsv")

A.sp_ALSP_WGD_paralogs_peak2 <- A.sp_ALSP %>% 
  filter(Ks > 0.772909) %>% 
  filter(Ks < 1.930947) %>% 
  select(X)

write.table(A.sp_ALSP_WGD_paralogs_peak2, quote = F, file = "Alsophila_sp_ALSP_WGD_paralogs_peak2.tsv")

#### Alsophila spinulosa ASLI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.spinulosa_ALSI <- read.delim("ks_distributions/Alsophila_spinulosa.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.spinulosa_ALSI_filt <- A.spinulosa_ALSI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.spinulosa_ALSI_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.spinulosa_ALSI_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Alsophila spinulosa (ALSI)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.spinulosa_ALSI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.spinulosa_ALSI_filt.0 <- A.spinulosa_ALSI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.spinulosa_ALSI_normalmixEM <- normalmixEM(A.spinulosa_ALSI_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3 
summary(A.spinulosa_ALSI_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.365029 0.408358 0.226614
#mu     0.278668 1.293813 2.913275
#sigma  0.146794 0.557326 0.615491
#loglik at estimate:  -12402.86 

# Generate list of paralogs 

A.spinulosa_ALSI_WGD_paralogs_peak1 <- A.spinulosa_ALSI %>% 
  filter(Ks > 0.13187) %>% 
  filter(Ks < 0.425462) %>% 
  select(X)

write.table(A.spinulosa_ALSI_WGD_paralogs_peak1, quote = F, file = "Alsophila_spinulosa_ALSI_WGD_paralogs_peak1.tsv")

A.spinulosa_ALSI_WGD_paralogs_peak2 <- A.spinulosa_ALSI %>% 
  filter(Ks > 0.736487) %>% 
  filter(Ks < 1.851139) %>% 
  select(X)

write.table(A.spinulosa_ALSI_WGD_paralogs_peak2, quote = F, file = "Alsophila_spinulosa_ALSI_WGD_paralogs_peak2.tsv")

#### Anemia phyllitidis ANPH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.phyllitidis_ANPH <- read.delim("ks_distributions/Anemia_phyllitidis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.phyllitidis_ANPH_filt <- A.phyllitidis_ANPH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.phyllitidis_ANPH_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.phyllitidis_ANPH_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Anemia phyllitidis (ANPH)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.phyllitidis_ANPH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.phyllitidis_ANPH_filt.0 <- A.phyllitidis_ANPH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.phyllitidis_ANPH_normalmixEM <- normalmixEM(A.phyllitidis_ANPH_filt.0$Ks, k=2)
# Fits to 2 WGD peaks, ignore comp 3 
summary(A.phyllitidis_ANPH_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1677206 0.832279
#mu     0.0864382 2.057770
#sigma  0.0764292 0.987752
#loglik at estimate:  -7345.379 

# Generate list of paralogs 

A.phyllitidis_ANPH_WGD_paralogs_peak1 <- A.phyllitidis_ANPH %>% 
  filter(Ks > 1.070018) %>% 
  filter(Ks < 3.045522) %>% 
  select(X)

write.table(A.phyllitidis_ANPH_WGD_paralogs_peak1, quote = F, file = "Anemia_phyllitidis_ANPH_WGD_paralogs_peak1.tsv")

#### Anemia tomenatosa CQPW #####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.tomentosa_CQPW <- read.delim("ks_distributions/Anemia_tomentosa_CQPW.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.tomentosa_CQPW_filt <- A.tomentosa_CQPW %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.tomentosa_CQPW_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.tomentosa_CQPW_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Anemia tomentosa (CQPW)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.tomentosa_CQPW.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.tomentosa_CQPW_filt.0 <- A.tomentosa_CQPW_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.tomentosa_CQPW_normalmixEM <- normalmixEM(A.tomentosa_CQPW_filt.0$Ks)
summary(A.tomentosa_CQPW_normalmixEM)

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.1877785 0.812222
#mu     0.0358278 1.908926
#sigma  0.0322500 1.064283
#loglik at estimate:  -7172.595

# Generate list of paralogs 

A.tomentosa_CQPW_WGD_paralogs_peak1 <- A.tomentosa_CQPW %>% 
  filter(Ks > 0.844643) %>% 
  filter(Ks < 2.973209) %>% 
  select(X)

write.table(A.tomentosa_CQPW_WGD_paralogs_peak1, quote = F, file = "Anemia_tomentosa_CQPW_WGD_paralogs_peak1.tsv")

#### Angiopteris fokiensis ANFK ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.fokiensis_ANFK <- read.delim("ks_distributions/Angiopteris_fokiensis_ANFK.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.fokiensis_ANFK_filt <- A.fokiensis_ANFK %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.fokiensis_ANFK_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.fokiensis_ANFK_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Angiopteris fokiensis (ANFK)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.fokiensis_ANFK.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.fokiensis_ANFK_filt.0 <- A.fokiensis_ANFK_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.fokiensis_ANFK_CQPW_normalmixEM <- normalmixEM(A.fokiensis_ANFK_filt.0$Ks, k=2)
#only fits to one WGD, ignore second comp
summary(A.fokiensis_ANFK_CQPW_normalmixEM)

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.538842 0.461158
#mu     0.784968 2.466250
#sigma  0.404386 0.820076
#loglik at estimate:  -7331.865

# Generate list of paralogs 

A.fokiensis_ANFK_WGD_paralogs_peak1 <- A.fokiensis_ANFK %>% 
  filter(Ks > 0.380582) %>% 
  filter(Ks < 1.189354) %>% 
  select(X)

write.table(A.fokiensis_ANFK_WGD_paralogs_peak1, quote = F, file = "Angiopteris_fokiensis_ANFK_WGD_paralogs_peak1.tsv")

#### Angiopteris fokiensis ANFO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.fokiensis_ANFO <- read.delim("ks_distributions/Angiopteris_fokiensis_ANFO.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.fokiensis_ANFO_filt <- A.fokiensis_ANFO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.fokiensis_ANFO_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.fokiensis_ANFO_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Angiopteris fokiensis (ANFO)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.fokiensis_ANFO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.fokiensis_ANFO_filt.0 <- A.fokiensis_ANFO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.fokiensis_ANFO_normalmixEM <- normalmixEM(A.fokiensis_ANFO_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 1 
summary(A.fokiensis_ANFO_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.434825 0.565175
#mu     2.459517 0.768211
#sigma  0.848342 0.425526
#loglik at estimate:  -8763.795  

# Generate list of paralogs 

A.fokiensis_ANFO_WGD_paralogs <- A.fokiensis_ANFO %>% 
  filter(Ks > 0.342685) %>% 
  filter(Ks < 1.193737) %>% 
  select(X)

write.table(A.fokiensis_ANFO_WGD_paralogs, quote = F, file = "Angiopteris_fokiensis_ANFO_WGD_paralogs.tsv")

#### Antrophym callifolium ANCA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.callifolium_ANCA <- read.delim("ks_distributions/Antrophyum_callifolium.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.callifolium_ANCA_filt <- A.callifolium_ANCA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.callifolium_ANCA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.callifolium_ANCA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Antrophym callifolium (ANCA)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.callifolium_ANCA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.callifolium_ANCA_filt.0 <- A.callifolium_ANCA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.callifolium_ANCA_normalmixEM <- normalmixEM(A.callifolium_ANCA_filt.0$Ks, k=2)
#only fits to one WGD, ignore first comp
summary(A.callifolium_ANCA_normalmixEM)

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.497354 0.502646
#mu     2.466864 0.530737
#sigma  0.872846 0.303595
#loglik at estimate:  -4811.377 

# Generate list of paralogs 

A.callifolium_ANCA_WGD_paralogs <- A.callifolium_ANCA %>% 
  filter(Ks > 0.227142) %>% 
  filter(Ks < 0.834332) %>% 
  select(X)

write.table(A.callifolium_ANCA_WGD_paralogs, quote = F, file = "Antrophym_callifolium_ANCA_WGD_paralogs.tsv")


#### Arachniodes nigrospinosa ARNI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.nigrospinosa_ARNI <- read.delim("ks_distributions/Arachniodes_nigrospinosa.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.nigrospinosa_ARNI_filt <- A.nigrospinosa_ARNI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.nigrospinosa_ARNI_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.nigrospinosa_ARNI_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Arachniodes nigrospinosa (ARNI)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.nigrospinosa_ARNI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.nigrospinosa_ARNI_filt.0 <- A.nigrospinosa_ARNI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.nigrospinosa_ARNI_normalmixEM <- normalmixEM(A.nigrospinosa_ARNI_filt.0$Ks, k=2)
#only fits to one WGD, ignore first comp
summary(A.nigrospinosa_ARNI_normalmixEM)

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.121206 0.878794
#mu     0.198461 2.198655
#sigma  0.153162 0.944954
#loglik at estimate:  -7651.55 

# Generate list of paralogs 

A.nigrospinosa_ARNI_WGD_paralogs <- A.nigrospinosa_ARNI %>% 
  filter(Ks > 1.253701) %>% 
  filter(Ks < 3.143609) %>% 
  select(X)

write.table(A.nigrospinosa_ARNI_WGD_paralogs, quote = F, file = "Arachniodes_nigrospinosa_ARNI_WGD_paralogs.tsv")

#### Argyrochosma nivea XDDT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.nivea_XDDT <- read.delim("ks_distributions/Argyrochosma_nivea_XDDT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.nivea_XDDT_filt <- A.nivea_XDDT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.nivea_XDDT_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.nivea_XDDT_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Argyrochosma nivea (XDDT)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.nivea_XDDT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.nivea_XDDT_filt.0 <- A.nivea_XDDT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.nivea_XDDT_BMJR_normalmixEM <- normalmixEM(A.nivea_XDDT_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.raddianum_BMJR_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.531089 0.468911
#mu     0.286095 2.368023
#sigma  0.308078 0.911571
#loglik at estimate:  -10397.71 

# Generate list of paralogs 

A.nivea_XDDT_WGD_paralogs <- A.nivea_XDDT %>% 
  filter(Ks > 1.456452) %>% 
  filter(Ks < 3.279594) %>% 
  select(X)

write.table(A.nivea_XDDT_WGD_paralogs, quote = F, file = "Argyrochosma_nivea_XDDT_WGD_paralogs.tsv")

#### Arthromeris lungtanuensis ARLU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.luntanuensis_ARLU <- read.delim("ks_distributions/Arthromeris_lungtauensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.luntanuensis_ARLU_filt <- A.luntanuensis_ARLU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.luntanuensis_ARLU_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.luntanuensis_ARLU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Arthromeris lungtauensis (ARLU)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.lungtauensis_ARLU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.luntanuensis_ARLU_filt.0 <- A.luntanuensis_ARLU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.luntanuensis_ARLU_normalmixEM <- normalmixEM(A.luntanuensis_ARLU_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 1 
summary(A.luntanuensis_ARLU_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.136609 0.863391
#mu     0.219441 2.254657
#sigma  0.165792 0.943773
#loglik at estimate:  -9042.911 

# Generate list of paralogs 

A.luntanuensis_ARLU_WGD_paralogs <- A.luntanuensis_ARLU %>% 
  filter(Ks > 1.310884) %>% 
  filter(Ks < 3.19843) %>% 
  select(X)

write.table(A.luntanuensis_ARLU_WGD_paralogs, quote = F, file = "Arthromeris_lungtauensis_ARLU_WGD_paralogs.tsv")

#### Arthropteris palisotii ARPI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.palisotii_ARPI<- read.delim("ks_distributions/Arthropteris_palisotii_ARPI.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.palisotii_ARPI_filt <- A.palisotii_ARPI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.palisotii_ARPI_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.palisotii_ARPI_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Arthropteris palisotii (ARPI)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.A.palisotii_ARPI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.palisotii_ARPI_filt.0 <- A.palisotii_ARPI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.palisotii_ARPI_normalmixEM <- normalmixEM(A.palisotii_ARPI_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.palisotii_ARPI_normalmixEM )
 
#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.162050 0.837950
#mu     0.128435 2.135616
#sigma  0.098917 0.982445
#loglik at estimate:  -7901.441

# Generate list of paralogs 

A.palisotii_ARPI_WGD_paralogs <- A.palisotii_ARPI %>% 
  filter(Ks > 1.153171) %>% 
  filter(Ks < 3.118061) %>% 
  select(X)

write.table(A.palisotii_ARPI_WGD_paralogs, quote = F, file = "Arthropteris_palisotii_ARPI_WGD_paralogs.tsv")

#### Arthropteris palisotii ARPA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.palisotii_ARPA<- read.delim("ks_distributions/Arthropteris_palisotii_ARPA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.palisotii_ARPA_filt <- A.palisotii_ARPA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.palisotii_ARPA_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.palisotii_ARPA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Arthropteris palisotii (ARPA)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.palisotii_ARPA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.palisotii_ARPA_filt.0 <- A.palisotii_ARPA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.palisotii_ARPA_normalmixEM <- normalmixEM(A.palisotii_ARPA_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.palisotii_ARPA_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1190146 0.880985
#mu     0.0903346 2.183167
#sigma  0.0816364 1.002653
#loglik at estimate:  -6323.15 

# Generate list of paralogs 

A.palisotii_ARPA_WGD_paralogs <- A.palisotii_ARPA %>% 
  filter(Ks > 1.180514) %>% 
  filter(Ks < 3.18582) %>% 
  select(X)

write.table(A.palisotii_ARPA_WGD_paralogs, quote = F, file = "Arthropteris_palisotii_ARPA_WGD_paralogs.tsv")

#### Arthropteris repens ARRE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.repens_ARRE<- read.delim("ks_distributions/Arthropteris_repens.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.repens_ARRE_filt <- A.repens_ARRE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.repens_ARRE_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.repens_ARRE_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Arthropteris repens (ARRE)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.repens_ARRE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.repens_ARRE_filt.0 <- A.repens_ARRE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.repens_ARRE_normalmixEM <- normalmixEM(A.repens_ARRE_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.repens_ARRE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.115269 0.884731
#mu     0.164118 2.208426
#sigma  0.119196 0.955563
#loglik at estimate:  -8238.681

# Generate list of paralogs 

A.repens_ARRE_WGD_paralogs <- A.repens_ARRE %>% 
  filter(Ks > 1.252863) %>% 
  filter(Ks < 3.163989) %>% 
  select(X)

write.table(A.repens_ARRE_WGD_paralogs, quote = F, file = "Arthropteris_repens_ARRE_WGD_paralogs.tsv")

#### Asplenium loriceum ASLO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.loriceum_ASLO<- read.delim("ks_distributions/Asplenium_loriceum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.loriceum_ASLO_filt <- A.loriceum_ASLO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.loriceum_ASLO_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.loriceum_ASLO_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Asplenium_loriceum (ASLO)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.loriceum_ASLO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.loriceum_ASLO_filt.0 <- A.loriceum_ASLO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.loriceum_ASLO_normalmixEM <- normalmixEM(A.loriceum_ASLO_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore third comp 
summary(A.loriceum_ASLO_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2    comp 3
#lambda 0.486552 0.438621 0.0748268
#mu     0.343965 2.106608 3.5844725
#sigma  0.182735 0.804226 0.2585708
#loglik at estimate:  -8265.614 

# Generate list of paralogs 

A.loriceum_ASLO_WGD_paralogs_peak1 <- A.loriceum_ASLO %>% 
  filter(Ks > 0.16123) %>% 
  filter(Ks < 0.5267) %>% 
  select(X)

write.table(A.loriceum_ASLO_WGD_paralogs_peak1, quote = F, file = "Asplenium_loriceum_ASLO_WGD_paralogs_peak1.tsv")

A.loriceum_ASLO_WGD_paralogs_peak2 <- A.loriceum_ASLO %>% 
  filter(Ks > 1.302382) %>% 
  filter(Ks < 2.910834) %>% 
  select(X)

write.table(A.loriceum_ASLO_WGD_paralogs_peak2, quote = F, file = "Asplenium_loriceum_ASLO_WGD_paralogs_peak2.tsv")

#### Asplenium pekinense ASPE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.pekinense_ASPE <- read.delim("ks_distributions/Asplenium_pekinense.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.pekinense_ASPE_filt <- A.pekinense_ASPE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.pekinense_ASPE_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.pekinense_ASPE_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Asplenium pekinense (ASPE)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.pekinense_ASPE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.pekinense_ASPE_filt.0 <- A.pekinense_ASPE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.pekinense_ASPE_normalmixEM <- normalmixEM(A.pekinense_ASPE_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(A.pekinense_ASPE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.200523 0.799477
#mu     0.262895 2.277884
#sigma  0.219586 0.948738
#loglik at estimate:  -9192.656  

#Generate list of paralogs 

A.pekinense_ASPE_WGD_paralogs <- A.pekinense_ASPE %>% 
  filter(Ks > 1.329146) %>% 
  filter(Ks < 3.226622) %>% 
  select(X)

write.table(A.pekinense_ASPE_WGD_paralogs, quote = F, file = "Asplenium_pekinense_ASPE_WGD_paralogs.tsv")

#### Asplenium platyneuron KJZG  #####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.platyneuron_KJZG <- read.delim("ks_distributions/Asplenium_platyneuron_KJZG.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.platyneuron_KJZG_filt <- A.platyneuron_KJZG %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.platyneuron_KJZG_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.platyneuron_KJZG_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Asplenium platyneuron (KJZG)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.platyneuron_KJZG.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.platyneuron_KJZG_filt.0 <- A.platyneuron_KJZG_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.platyneuron_KJZG_normalmixEM <- normalmixEM(A.platyneuron_KJZG_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.platyneuron_KJZG_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.132087 0.867913
#mu     0.138961 2.243912
#sigma  0.109485 0.967108
#loglik at estimate:  -6442.227

# Generate list of paralogs 

A.platyneuron_KJZG_WGD_paralogs <- A.platyneuron_KJZG %>% 
  filter(Ks > 1.276804) %>% 
  filter(Ks < 3.21102) %>% 
  select(X)

write.table(A.platyneuron_KJZG_WGD_paralogs, quote = F, file = "Asplenium_platyneuron_KJZG_WGD_paralogs.tsv")

#### Asplenium polyodon ASPY ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.polyodon_ASPY <- read.delim("ks_distributions/Asplenium_polyodon.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.polyodon_ASPY_filt <- A.polyodon_ASPY %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.polyodon_ASPY_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.polyodon_ASPY_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Asplenium polyodon (ASPY)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.polyodon_ASPY.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.polyodon_ASPY_filt.0 <- A.polyodon_ASPY_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.polyodon_ASPY_normalmixEM <- normalmixEM(A.polyodon_ASPY_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore third comp 
summary(A.polyodon_ASPY_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.430876 0.419815 0.149309
#mu     1.939654 0.407834 3.406605
#sigma  0.741696 0.222079 0.360394
#loglik at estimate:  -6565.283 

# Generate list of paralogs 

A.polyodon_ASPY_WGD_paralogs_peak1 <- A.polyodon_ASPY %>% 
  filter(Ks > 0.185755) %>% 
  filter(Ks < 0.629913) %>% 
  select(X)

write.table(A.polyodon_ASPY_WGD_paralogs_peak1, quote = F, file = "Asplenium_polyodon_ASPY_WGD_paralogs_peak1.tsv")

A.polyodon_ASPY_WGD_paralogs_peak2 <- A.polyodon_ASPY %>% 
  filter(Ks > 1.197958) %>% 
  filter(Ks < 2.68135) %>% 
  select(X)

write.table(A.polyodon_ASPY_WGD_paralogs_peak2, quote = F, file = "Asplenium_polyodon_ASPY_WGD_paralogs_peak2.tsv")

#### Asplenium ruprechtii ASRU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.ruprechtii_ASRU <- read.delim("ks_distributions/Asplenium_ruprechtii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.ruprechtii_ASRU_filt <- A.ruprechtii_ASRU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.ruprechtii_ASRU_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.ruprechtii_ASRU_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Asplenium ruprechtii (ASRU)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.ruprechtii_ASRU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.ruprechtii_ASRU_filt.0 <- A.ruprechtii_ASRU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.ruprechtii_ASRU_normalmixEM <- normalmixEM(A.ruprechtii_ASRU_filt.0$Ks, k=2)
# Fits to 2 WGD peaks, ignore third comp 
summary(A.ruprechtii_ASRU_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.124384 0.875616
#mu     0.210711 2.341796
#sigma  0.171575 0.919470
#loglik at estimate:  -5997.096

# Generate list of paralogs 

A.ruprechtii_ASRU_WGD_paralogs <- A.ruprechtii_ASRU %>% 
  filter(Ks > 1.422326) %>% 
  filter(Ks < 3.261266) %>% 
  select(X)

write.table(A.ruprechtii_ASRU_WGD_paralogs, quote = F, file = "Asplenium_ruprechtii_ASRU_WGD_paralogs.tsv")

#### Athyrium decurrenti-alatum ATDA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.decurrenti_ATDA <- read.delim("ks_distributions/Athyrium_decurrenti-alatum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.decurrenti_ATDA_filt <- A.decurrenti_ATDA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.decurrenti_ATDA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.decurrenti_ATDA_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Athyrium decurrenti-alatum (ATDA)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.decurrenti_ATDA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.decurrenti_ATDA_filt.0 <- A.decurrenti_ATDA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.decurrenti_ATDA_normalmixEM <- normalmixEM(A.decurrenti_ATDA_filt.0$Ks, k=2)
# Fits to 1 WGD peaks, ignore third comp 
summary(A.decurrenti_ATDA_normalmixEM )

#summary of normalmixEM object:
#           comp 1   comp 2
#lambda 0.0776308 0.922369
#mu     0.1241702 2.022094
#sigma  0.1009716 0.992390
#loglik at estimate:  -5943.68

# Generate list of paralogs 

A.decurrenti_ATDA_WGD_paralogs <- A.decurrenti_ATDA %>% 
  filter(Ks > 1.029704) %>% 
  filter(Ks < 3.014484) %>% 
  select(X)

write.table(A.decurrenti_ATDA_WGD_paralogs, quote = F, file = "Athyrium_decurrenti-alatum_ADTA_WGD_paralogs.tsv")

#### Athyrium filix-femina URCP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.filix_femina_URCP <- read.delim("ks_distributions/Athyrium_filix-femina_URCP.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.filix_femina_URCP_filt <- A.filix_femina_URCP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.filix_femina_URCP_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.filix_femina_URCP_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Athyrium filix-femina (URCP)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.filix-femina_URCP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.filix_femina_URCP_filt.0 <- A.filix_femina_URCP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.filix_femina_URCP_normalmixEM <- normalmixEM(A.filix_femina_URCP_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(A.filix_femina_URCP_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.1111353 0.888865
#mu     0.1121738 2.105979
#sigma  0.0875735 0.985435
#loglik at estimate:  -8858.296

# Generate list of paralogs 

A.filix_femina_URCP_WGD_paralogs <- A.filix_femina_URCP %>% 
  filter(Ks > 1.120544) %>% 
  filter(Ks < 3.091414) %>% 
  select(X)

write.table(A.filix_femina_URCP_WGD_paralogs, quote = F, file = "Athyrium_filix_femina_URCP_WGD_paralogs.tsv")

#### Athyrium iseanum ATIS ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.isaenum_ATIS <- read.delim("ks_distributions/Athyrium_iseanum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.isaenum_ATIS_filt <- A.isaenum_ATIS %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.isaenum_ATIS_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.isaenum_ATIS_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Athyrium isaenum (ATIS)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.isaenum_ATIS.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.isaenum_ATIS_filt.0 <- A.isaenum_ATIS_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.isaenum_ATIS_normalmixEM <- normalmixEM(A.isaenum_ATIS_filt.0$Ks, k=2)
# Fits to 1 WGD peaks, ignore third comp 
summary(A.isaenum_ATIS_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1150291 0.884971
#mu     0.0951026 2.041746
#sigma  0.0817642 1.001470
#loglik at estimate:  -8022.275

# Generate list of paralogs 

A.isaenum_ATIS_WGD_paralogs <- A.isaenum_ATIS %>% 
  filter(Ks > 1.040276) %>% 
  filter(Ks < 3.043216) %>% 
  select(X)

write.table(A.isaenum_ATIS_WGD_paralogs, quote = F, file = "Athyrium_isaenum_ATIS_WGD_paralogs.tsv")

#### Azolla cf. caroliniana CVEG ####

A.cf.caroliniana_CVEG <- read.delim("ks_distributions/Azolla_caroliniana_CVEG.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

A.cf.caroliniana_CVEG_filt <- A.cf.caroliniana_CVEG %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.cf.caroliniana_CVEG_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.cf.caroliniana_CVEG_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Azolla cf. caroliniana (CVEG)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.cf.caroliniana_CVEG.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.cf.caroliniana_CVEG_filt.0 <- A.cf.caroliniana_CVEG_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

#No WGD peak in Ks plot

#### Azolla pinnata AZPN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.pinnata_AZPN <- read.delim("ks_distributions/Azolla_pinnata.cds.ks.tsv")

# Filter out saturated duplicates (Ks < 4)

A.pinnata_AZPN_filt <- A.pinnata_AZPN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.pinnata_AZPN_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.pinnata_AZPN_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Azolla pinnata (AZPN)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.pinnata_AZPN.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Azolla pinnata AZPI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

A.pinnata_AZPI <- read.delim("ks_distributions/Azolla_pinnata_AZPI.cds.ks.tsv")

# Filter out saturated duplicates (Ks < 4)

A.pinnata_AZPI_filt <- A.pinnata_AZPI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(A.pinnata_AZPI_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(A.pinnata_AZPI_filt, mapping = aes(x=Ks, ..scaled..*3000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Azolla pinnata (AZPI)") + theme(plot.title = element_text(face = "italic"))

ggsave("A.pinnata_AZPI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

A.pinnata_AZPI_filt.0 <- A.pinnata_AZPI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

A.pinnata_AZPI_normalmixEM <- normalmixEM(A.pinnata_AZPI_filt.0$Ks, k=3)
# Fits to 1 WGD peak, ignore third comp 
summary(A.pinnata_AZPI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.2937287 0.223440 0.482831
#mu     0.1061189 0.774378 2.481781
#sigma  0.0837476 0.378093 0.833703
#loglik at estimate:  -6341.191 

# Generate list of paralogs 

A.pinnata_AZPI_WGD_paralogs <- A.pinnata_AZPI %>% 
  filter(Ks > 0.396285) %>% 
  filter(Ks < 1.15471) %>% 
  select(X)

write.table(A.pinnata_AZPI_WGD_paralogs, quote = F, file = "Azolla_pinnata_AZPI_WGD_paralogs.tsv")

#### Blechnopsis orientalis BLOR ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

B.orientalis_BLOR <- read.delim("ks_distributions/Blechnopsis_orientalis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

B.orientalis_BLOR_filt <- B.orientalis_BLOR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(B.orientalis_BLOR_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(B.orientalis_BLOR_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Blechnopsis orientalis (BLOR)") + theme(plot.title = element_text(face = "italic"))

ggsave("B.orientalis_BLOR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

B.orientalis_BLOR_filt.0 <- B.orientalis_BLOR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

B.orientalis_BLOR_normalmixEM <- normalmixEM(B.orientalis_BLOR_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(B.orientalis_BLOR_normalmixEM )

#summary of normalmixEM object:
#          comp 1   comp 2
#lambda 0.114361 0.885639
#mu     0.171810 2.133826
#sigma  0.140442 0.971520
#loglik at estimate:  -8134.456 

# Generate list of paralogs 

B.orientalis_BLOR_WGD_paralogs <- B.orientalis_BLOR %>% 
  filter(Ks > 1.162306) %>% 
  filter(Ks < 3.105346) %>% 
  select(X)

write.table(B.orientalis_BLOR_WGD_paralogs, quote = F, file = "Blechnopsis_orientalis_BLOR_WGD_paralogs.tsv")

#### Blechnum spicant AFPO  ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

B.spicant_AFPO <- read.delim("ks_distributions/Blechnum_spicant_AFPO.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

B.spicant_AFPO_filt <- B.spicant_AFPO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(B.spicant_AFPO_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(B.spicant_AFPO_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Blechnum spicant (AFPO)") + theme(plot.title = element_text(face = "italic"))

ggsave("B.spicant_AFPO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

B.spicant_AFPO_filt.0 <- B.spicant_AFPO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

B.spicant_AFPO_normalmixEM <- normalmixEM(B.spicant_AFPO_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(B.spicant_AFPO_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.897903 0.1020967
#mu     2.032521 0.1364486
#sigma  0.996930 0.0987764
#loglik at estimate:  -9238.625

# Generate list of paralogs 

B.spicant_AFPO_WGD_paralogs <- B.spicant_AFPO %>% 
  filter(Ks > 1.035591) %>% 
  filter(Ks < 3.029451) %>% 
  select(X)

write.table(B.spicant_AFPO_WGD_paralogs, quote = F, file = "Blechnum_spicant_AFPO_WGD_paralogs.tsv")

#### Bolbitis appendiculata BOAP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

B.appendiculata_BOAP <- read.delim("ks_distributions/Bolbitis_appendiculata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

B.appendiculata_BOAP_filt <- B.appendiculata_BOAP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(B.appendiculata_BOAP_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(B.appendiculata_BOAP_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Bolbitis appendiculata (BOAP)") + theme(plot.title = element_text(face = "italic"))

ggsave("B.appendiculata_BOAP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

B.appendiculata_BOAP_filt.0 <- B.appendiculata_BOAP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

B.appendiculata_BOAP_normalmixEM <- normalmixEM(B.appendiculata_BOAP_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(B.appendiculata_BOAP_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.174508 0.825492
#mu     0.244495 2.245405
#sigma  0.191930 0.955766
#loglik at estimate:  -11966.93 

# Generate list of paralogs 

B.appendiculata_BOAP_WGD_paralogs <- B.appendiculata_BOAP %>% 
  filter(Ks > 1.289639) %>% 
  filter(Ks < 3.201171) %>% 
  select(X)

write.table(B.appendiculata_BOAP_WGD_paralogs, quote = F, file = "Bolbitis_appendiculata_BOAP_WGD_paralogs.tsv")

#### Bolbitis heteroclita BOHE ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

B.heteroclita_BOHE <- read.delim("ks_distributions/Bolbitis_heteroclita.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

B.heteroclita_BOHE_filt <- B.heteroclita_BOHE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(B.heteroclita_BOHE_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(B.heteroclita_BOHE_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Bolbitis heteroclita (BOHE)") + theme(plot.title = element_text(face = "italic"))

ggsave("B.heteroclita_BOHE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

B.heteroclita_BOHE_filt.0 <- B.heteroclita_BOHE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

B.heteroclita_BOHE_normalmixEM <- normalmixEM(B.heteroclita_BOHE_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(B.heteroclita_BOHE_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.167423 0.832577
#mu     0.238366 2.161716
#sigma  0.189460 0.964583
#loglik at estimate:  -8464.497 

# Generate list of paralogs 

B.heteroclita_BOHE_WGD_paralogs <- B.heteroclita_BOHE %>% 
  filter(Ks > 1.197133) %>% 
  filter(Ks < 3.126299) %>% 
  select(X)

write.table(B.heteroclita_BOHE_WGD_paralogs, quote = F, file = "Bolbitis_heteroclita_BOHE_WGD_paralogs.tsv")

#### Bosmania membranacea BOME ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

B.membranacea_BOME <- read.delim("ks_distributions/Bosmania_membranacea.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

B.membranacea_BOME_filt <- B.membranacea_BOME %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(B.membranacea_BOME_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(B.membranacea_BOME_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Bosmania membranacea (BOME)") + theme(plot.title = element_text(face = "italic"))

ggsave("B.membranacea_BOME.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

B.membranacea_BOME_filt.0 <- B.membranacea_BOME_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

B.membranacea_BOME_normalmixEM <- normalmixEM(B.membranacea_BOME_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(B.membranacea_BOME_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.118820 0.881180
#mu     0.224335 2.227340
#sigma  0.164540 0.942184
#loglik at estimate:  -7376.135 

# Generate list of paralogs 

B.membranacea_BOME_WGD_paralogs <- B.membranacea_BOME %>% 
  filter(Ks > 1.285156) %>% 
  filter(Ks < 3.169524) %>% 
  select(X)

write.table(B.membranacea_BOME_WGD_paralogs, quote = F, file = "Bosmania_membranacea_BOME_WGD_paralogs.tsv")

#### Callistopteris apifolia CAAP ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.apifolia_CAAP <- read.delim("ks_distributions/Callistopteris_apifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.apifolia_CAAP_filt <- C.apifolia_CAAP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.apifolia_CAAP_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.apifolia_CAAP_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Callistopteris apifolia (CAAP)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.apifolia_CAAP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.apifolia_CAAP_filt.0 <- C.apifolia_CAAP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.apifolia_CAAP_normalmixEM <- normalmixEM(C.apifolia_CAAP_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(C.apifolia_CAAP_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.463382 0.536618
#mu     0.817481 2.670585
#sigma  0.522586 0.739348
#loglik at estimate:  -6864.125

# Generate list of paralogs 

C.apifolia_CAAP_WGD_paralogs <- C.apifolia_CAAP %>% 
  filter(Ks > 0.294895) %>% 
  filter(Ks < 1.340067) %>% 
  select(X)

write.table(C.apifolia_CAAP_WGD_paralogs, quote = F, file = "Callistopteris_apifolia_CAAP_WGD_paralogs.tsv")

#### Cephalomanes javanicum CEJA ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.javanicum_CEJA <- read.delim("ks_distributions/Cephalomanes_javanicum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.javanicum_CEJA_filt <- C.javanicum_CEJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.javanicum_CEJA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.javanicum_CEJA_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cephalomanes javanicum (CEJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.javanicum_CEJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.javanicum_CEJA_filt.0 <- C.javanicum_CEJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.javanicum_CEJA_normalmixEM <- normalmixEM(C.javanicum_CEJA_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(C.javanicum_CEJA_normalmixEM )
#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.550505 0.449495
#mu     1.018244 2.872027
#sigma  0.637858 0.640736
#loglik at estimate:  -6104.876  

# Generate list of paralogs 

C.javanicum_CEJA_WGD_paralogs <- C.javanicum_CEJA %>% 
  filter(Ks > 0.380386) %>% 
  filter(Ks < 1.656102) %>% 
  select(X)

write.table(C.javanicum_CEJA_WGD_paralogs, quote = F, file = "Cephalomanes_javanicum_CEJA_WGD_paralogs.tsv")

#### Ceratopteris thalictroides CETH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.thalictroides_CETH <- read.delim("ks_distributions/Ceratopteris_thalicotroidesa.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.thalictroides_CETH_filt <- C.thalictroides_CETH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.thalictroides_CETH_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.thalictroides_CETH_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ceratopteris thalictroides (CETH)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.thalictroides_CETH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.thalictroides_CETH_filt.0 <- C.thalictroides_CETH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.thalictroides_CETH_normalmixEM <- normalmixEM(C.thalictroides_CETH_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(C.thalictroides_CETH_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.657165 0.342835
#mu     1.228955 3.040661
#sigma  0.738555 0.546495
#loglik at estimate:  -8263.438

# Generate list of paralogs 

C.thalictroides_CETH_WGD_paralogs <- C.thalictroides_CETH %>% 
  filter(Ks > 0.4904) %>% 
  filter(Ks < 1.96751) %>% 
  select(X)

write.table(C.thalictroides_CETH_WGD_paralogs, quote = F, file = "Ceratopteris_thalictroides_CETH_WGD_paralogs.tsv")

#### Cheilanthes chusana CHCH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.chusana_CHCH <- read.delim("ks_distributions/Cheilanthes_chusana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.chusana_CHCH_filt <- C.chusana_CHCH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.chusana_CHCH_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.chusana_CHCH_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cheilanthes chusana (CHCH)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.chusana_CHCH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.chusana_CHCH_filt.0 <- C.chusana_CHCH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.chusana_CHCH_normalmixEM <- normalmixEM(C.chusana_CHCH_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(C.chusana_CHCH_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.107148 0.892852
#mu     0.155422 2.357101
#sigma  0.128028 0.939615
#loglik at estimate:  -7369.894 

# Generate list of paralogs 

C.chusana_CHCH_WGD_paralogs <- C.chusana_CHCH %>% 
  filter(Ks > 1.417486) %>% 
  filter(Ks < 3.296716) %>% 
  select(X)

write.table(C.chusana_CHCH_WGD_paralogs, quote = F, file = "Cheilanthes_chusana_CHCH_WGD_paralogs.tsv")

#### Cheilanthes nitidula CHNI ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.nitidula_CHNI <- read.delim("ks_distributions/Cheilanthes_nitidula.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.nitidula_CHNI_filt <- C.nitidula_CHNI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.nitidula_CHNI_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.nitidula_CHNI_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cheilanthes nitidula (CHNI)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.nitidula_CHNI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.nitidula_CHNI_filt.0 <- C.nitidula_CHNI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.nitidula_CHNI_normalmixEM <- normalmixEM(C.nitidula_CHNI_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(C.nitidula_CHNI_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1254747 0.874525
#mu     0.1188589 2.227769
#sigma  0.0963706 0.964996
#loglik at estimate:  -5122.616   

# Generate list of paralogs 

C.nitidula_CHNI_WGD_paralogs <- C.nitidula_CHNI %>% 
  filter(Ks > 1.262773) %>% 
  filter(Ks < 3.192765) %>% 
  select(X)

write.table(C.nitidula_CHNI_WGD_paralogs, quote = F, file = "Cheilanthes_nitidula_CHNI_WGD_paralogs.tsv")

#### Cheiropleuria bicuspis CHBI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.bicuspis_CHBI <- read.delim("ks_distributions/Cheiropleuria_bicuspis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.bicuspis_CHBI_filt <- C.bicuspis_CHBI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.bicuspis_CHBI_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.bicuspis_CHBI_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cheiropleuria bicuspis (CHBI)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.bicuspis_CHBI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.bicuspis_CHBI_filt.0 <- C.bicuspis_CHBI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.bicuspis_CHBI_normalmixEM <- normalmixEM(C.bicuspis_CHBI_filt.0$Ks, k=4)
# Fits to 1 WGD peak, comp 1 fits recent dups, comp 2 fits WGD, comps 3 and 4 are overfit
summary(C.bicuspis_CHBI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3   comp 4
#lambda 0.1172438 0.328932 0.448516 0.105308
#mu     0.0666120 0.768656 2.207356 3.556882
#sigma  0.0612058 0.365776 0.687945 0.269542
#loglik at estimate:  -7743.662 

# Generate list of paralogs 

C.bicuspis_CHBI_WGD_paralogs <- C.bicuspis_CHBI %>% 
  filter(Ks > 0.40288) %>% 
  filter(Ks < 1.134432) %>% 
  select(X)

write.table(C.bicuspis_CHBI_WGD_paralogs, quote = F, file = "Cheiropleuria_bicuspis_CHBI_WGD_paralogs.tsv")

#### Cheiropleuria integrifolia CHIN ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.integrifolia_CHIN <- read.delim("ks_distributions/Cheiropleuria_integrifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.integrifolia_CHIN_filt <- C.integrifolia_CHIN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.integrifolia_CHIN_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.integrifolia_CHIN_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cheiropleuria integrifolia (CHIN)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.integrifolia_CHIN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.integrifolia_CHIN_filt.0 <- C.integrifolia_CHIN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.integrifolia_CHIN_normalmixEM <- normalmixEM(C.integrifolia_CHIN_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2 (?)
summary(C.integrifolia_CHIN_normalmixEM )
#summary of normalmixEM object:
#comp 1   comp 2
#lambda 0.471082 0.528918
#mu     0.877318 2.668317
#sigma  0.560366 0.716393
#loglik at estimate:  -4683.211

# Generate list of paralogs 

C.integrifolia_CHIN_WGD_paralogs <- C.integrifolia_CHIN %>% 
  filter(Ks > 0.316952) %>% 
  filter(Ks < 1.437684) %>% 
  select(X)

write.table(C.integrifolia_CHIN_WGD_paralogs, quote = F, file = "Cheiropleuria_integrifolia_CHIN_WGD_paralogs.tsv")

#### Christella acuminata CHAC ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.acuminata_CHAC <- read.delim("ks_distributions/Christella_acuminata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.acuminata_CHAC_filt <- C.acuminata_CHAC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.acuminata_CHAC_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.acuminata_CHAC_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Christella acuminata (CHAC)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.acuminata_CHAC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.acuminata_CHAC_filt.0 <- C.acuminata_CHAC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.acuminata_CHAC_normalmixEM <- normalmixEM(C.acuminata_CHAC_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2 (?)
summary(C.acuminata_CHAC_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.111418 0.888582
#mu     0.117173 2.095508
#sigma  0.105128 0.978783
#loglik at estimate:  -9711.051 

# Generate list of paralogs 

C.acuminata_CHAC_WGD_paralogs <- C.acuminata_CHAC %>% 
  filter(Ks > 1.116725) %>% 
  filter(Ks < 3.074291) %>% 
  select(X)

write.table(C.acuminata_CHAC_WGD_paralogs, quote = F, file = "Christella_acuminata_CHAC_WGD_paralogs.tsv")

#### Christensenia aesculifolia CHAE ####

#Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.aesculifolia_CHAE <- read.delim("ks_distributions/Christensenia_aescuilfolira.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.aesculifolia_CHAE_filt <- C.aesculifolia_CHAE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.aesculifolia_CHAE_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.aesculifolia_CHAE_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Christensenia aesculifolia (CHAE)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.aesculifolia_CHAE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.aesculifolia_CHAE_filt.0 <- C.aesculifolia_CHAE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.aesculifolia_CHAE_normalmixEM <- normalmixEM(C.aesculifolia_CHAE_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2 (?)
summary(C.aesculifolia_CHAE_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.581646 0.418354
#mu     0.855877 2.675551
#sigma  0.542716 0.744162
#loglik at estimate:  -7081.883

# Generate list of paralogs 

C.aesculifolia_CHAE_WGD_paralogs <- C.aesculifolia_CHAE %>% 
  filter(Ks > 0.313161) %>% 
  filter(Ks < 1.398593) %>% 
  select(X)

write.table(C.aesculifolia_CHAE_WGD_paralogs, quote = F, file = "Christensenia_aesculifolia_CHAE_WGD_paralogs.tsv")

#### Cibotium barometz CIBZ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.barometz_CIBZ <- read.delim("ks_distributions/Cibotium_barometz_CIBZ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.barometz_CIBZ_filt <- C.barometz_CIBZ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.barometz_CIBZ_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.barometz_CIBZ_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cibotium barometz (CIBZ)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.barometz_CIBZ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.barometz_CIBZ_filt.0 <- C.barometz_CIBZ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.barometz_CIBZ_normalmixEM <- normalmixEM(C.barometz_CIBZ_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, comp 3 is overfit 
summary(C.barometz_CIBZ_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.367951 0.347523 0.284525
#mu     0.266647 1.177100 2.763231
#sigma  0.140374 0.492021 0.665939
#loglik at estimate:  -8676.918 

# Generate list of paralogs 

C.barometz_CIBZ_WGD_paralogs_peak1 <- C.barometz_CIBZ %>% 
  filter(Ks > 0.126273) %>% 
  filter(Ks < 0.407021) %>% 
  select(X)

write.table(C.barometz_CIBZ_WGD_paralogs_peak, quote = F, file = "Cibotium_barometz_CIBZ_WGD_paralogs_peak1.tsv")

C.barometz_CIBZ_WGD_paralogs_peak2 <- C.barometz_CIBZ %>% 
  filter(Ks > 0.685079) %>% 
  filter(Ks < 1.669121) %>% 
  select(X)

write.table(C.barometz_CIBZ_WGD_paralogs_peak2, quote = F, file = "Cibotium_barometz_CIBZ_WGD_paralogs_peak2.tsv")

#### Cibotium barometz CIBA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.barometz_CIBA <- read.delim("ks_distributions/Cibotium_barometz_CIBA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.barometz_CIBA_filt <- C.barometz_CIBA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.barometz_CIBA_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.barometz_CIBA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cibotium barometz (CIBA)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.barometz_CIBA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.barometz_CIBA_filt.0 <- C.barometz_CIBA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.barometz_CIBA_normalmixEM <- normalmixEM(C.barometz_CIBA_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, comp 3 is overfit 
summary(C.barometz_CIBA_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.396775 0.404059 0.199166
#mu     0.276120 1.337246 2.992645
#sigma  0.151427 0.591750 0.574601
#loglik at estimate:  -11244.98 

# Generate list of paralogs 

C.barometz_CIBA_WGD_paralogs_peak1 <- C.barometz_CIBA %>% 
  filter(Ks > 0.124693) %>% 
  filter(Ks < 0.427547) %>% 
  select(X)

write.table(C.barometz_CIBA_WGD_paralogs_peak1, quote = F, file = "Cibotium_barometz_CIBA_WGD_paralogs_peak1.tsv")

C.barometz_CIBA_WGD_paralogs_peak2 <- C.barometz_CIBA %>% 
  filter(Ks > 0.745496) %>% 
  filter(Ks < 1.928996) %>% 
  select(X)

write.table(C.barometz_CIBA_WGD_paralogs_peak2, quote = F, file = "Cibotium_barometz_CIBA_WGD_paralogs_peak2.tsv")

#### Coniogramme japonica COJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.japonica_COJA <- read.delim("ks_distributions/Coniogramme_japonica.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.japonica_COJA_filt <- C.japonica_COJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.japonica_COJA_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.japonica_COJA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Coniogramme japonica (COJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.japonica_COJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.japonica_COJA_filt.0 <- C.japonica_COJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.japonica_COJA_normalmixEM <- normalmixEM(C.japonica_COJA_filt.0$Ks, k=2)
# Fits to 1 WGD peaks
summary(C.japonica_COJA_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.0792384 0.920762
#mu     0.1035772 2.226581
#sigma  0.0886924 0.992731
#loglik at estimate:  -6449.452

# Generate list of paralogs 

C.japonica_COJA_WGD_paralogs <- C.japonica_COJA %>% 
  filter(Ks > 1.23385) %>% 
  filter(Ks < 3.219312) %>% 
  select(X)

write.table(C.japonica_COJA_WGD_paralogs, quote = F, file = "Coniogramme_japonica_COJA_WGD_paralogs.tsv")

#### Crepidomanes minutum CRMI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.minutum_CRMI <- read.delim("ks_distributions/Crepidomanes_minutum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.minutum_CRMI_filt <- C.minutum_CRMI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.minutum_CRMI_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.minutum_CRMI_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Crepidomanes minutum (CRMI)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.minutum_CRMI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.minutum_CRMI_filt.0 <- C.minutum_CRMI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.minutum_CRMI_normalmixEM <- normalmixEM(C.minutum_CRMI_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peaks, comp 2 is overfit? 
summary(C.minutum_CRMI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.560669 0.439331
#mu     1.100952 2.907270
#sigma  0.735041 0.621610
#loglik at estimate:  -4927.17  

# Generate list of paralogs 

C.minutum_CRMI_WGD_paralogs <- C.minutum_CRMI %>% 
  filter(Ks > 0.356911) %>% 
  filter(Ks < 1.835993) %>% 
  select(X)

write.table(C.minutum_CRMI_WGD_paralogs, quote = F, file = "Crepidomanes_minutum_CRMI_WGD_paralogs.tsv")

#### Cryptogramma acrostichoides WQML ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.acrostichoides_WQML <- read.delim("ks_distributions/Cryptogramma_acrostichoides_WQML.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.acrostichoides_WQML_filt <- C.acrostichoides_WQML %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.acrostichoides_WQML_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.acrostichoides_WQML_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cryptogramma acrostichoides (WQML)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.acrostichoides_WQML.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.acrostichoides_WQML_filt.0 <- C.acrostichoides_WQML_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.acrostichoides_WQML_normalmixEM <- normalmixEM(C.acrostichoides_WQML_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(C.acrostichoides_WQML_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.1125379 0.887462
#mu     0.1180393 2.261447
#sigma  0.0981544 0.974351
#loglik at estimate:  -8199.052

# Generate list of paralogs 

C.acrostichoides_WQML_WGD_paralogs <- C.acrostichoides_WQML %>% 
  filter(Ks > 1.287096) %>% 
  filter(Ks < 3.235798) %>% 
  select(X)

write.table(C.acrostichoides_WQML_WGD_paralogs, quote = F, file = "Cryptogramma_acrostichoides_WQML_WGD_paralogs.tsv")

#### Ctenitis subglandulosa CTSU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.subglandulosa_CTSU <- read.delim("ks_distributions/Ctenitis_subglanduolsa.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.subglandulosa_CTSU_filt <- C.subglandulosa_CTSU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.subglandulosa_CTSU_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.subglandulosa_CTSU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ctenitis subglandulosa (CTSU)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.subglandulosa_CTSU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.subglandulosa_CTSU_filt.0 <- C.subglandulosa_CTSU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.subglandulosa_CTSU_normalmixEM <- normalmixEM(C.subglandulosa_CTSU_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peaks, comp 2 is overfit? 
summary(C.subglandulosa_CTSU_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.0724638 0.927536
#mu     0.1139235 2.127755
#sigma  0.0942648 0.980042
#loglik at estimate:  -7982.992 

# Generate list of paralogs 

C.subglandulosa_CTSU_WGD_paralogs <- C.subglandulosa_CTSU %>% 
  filter(Ks > 1.147713) %>% 
  filter(Ks < 3.107797) %>% 
  select(X)

write.table(C.subglandulosa_CTSU_WGD_paralogs, quote = F, file = "Ctenitis_subglandulosa_CTSU_WGD_paralogs.tsv")

#### Culcita macrocarpa PNZO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.macrocarpa_PNZO <- read.delim("ks_distributions/Culcita_macrocarpa_PNZO.cds.ks.tsv")

# Filter out saturated duplicates (Ks < 4)

C.macrocarpa_PNZO_filt <- C.macrocarpa_PNZO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.macrocarpa_PNZO_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.macrocarpa_PNZO_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Culcita macrocarpa (PNZO)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.macrocarpa_PNZO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.macrocarpa_PNZO_filt.0 <- C.macrocarpa_PNZO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

# No evidence of WGD 

#### Cyathea spinulosa GANB #### 

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.spinulosa_GANB <- read.delim("ks_distributions/Cyathea_spinulosa_GANB.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.spinulosa_GANB_filt <- C.spinulosa_GANB %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.spinulosa_GANB_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.spinulosa_GANB_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cyathea (Aslophila) spinulosa (GANB)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.spinulosa_GANB.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.spinulosa_GANB_filt.0 <- C.spinulosa_GANB_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.spinulosa_GANB_normalmixEM <- normalmixEM(C.spinulosa_GANB_filt.0$Ks, k=3)
# Fits to 2 WGD peaks  
summary(C.spinulosa_GANB_normalmixEM )

# summary of normalmixEM object:
#          comp 1   comp 2   comp 3
#lambda 0.0524735 0.318910 0.628616
#mu     0.0361520 0.314880 1.906493
#sigma  0.0272385 0.141553 0.993765
#loglik at estimate:  -6682.747 

# Generate list of paralogs 

C.spinuolsa_GANB_WGD_paralogs_peak1 <- C.spinulosa_GANB %>% 
  filter(Ks > 0.173327) %>% 
  filter(Ks < 0.456433) %>% 
  select(X)

write.table(C.spinuolsa_GANB_WGD_paralogs_peak1, quote = F, file = "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak1.tsv")

C.spinuolsa_GANB_WGD_paralogs_peak2 <- C.spinulosa_GANB %>% 
  filter(Ks > 0.912728) %>% 
  filter(Ks < 2.900258) %>% 
  select(X)

write.table(C.spinuolsa_GANB_WGD_paralogs_peak2, quote = F, file = "Cyathea_Alsophila_spinuolsa_GANB_WGD_paralogs_peak2.tsv")

#### Cyclopeltis crenata CYCR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.crenata_CYCR <- read.delim("ks_distributions/Cyclopeltis_crenata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.crenata_CYCR_filt <- C.crenata_CYCR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.crenata_CYCR_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.crenata_CYCR_filt, mapping = aes(x=Ks, ..scaled..*400), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cyclopeltis crenata (CYCR)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.crenata_CYCR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.crenata_CYCR_filt.0 <- C.crenata_CYCR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.crenata_CYCR_normalmixEM <- normalmixEM(C.crenata_CYCR_filt.0$Ks, k=2, maxit = 2000)
# Fits to 2 WGD peaks, comp 3 is overfit 
summary(C.crenata_CYCR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.102396 0.897604
#mu     0.213476 2.200818
#sigma  0.170798 0.958479
#loglik at estimate:  -5290.492 

# Generate list of paralogs 

C.crenata_CYCR_WGD_paralogs_peak1 <- C.crenata_CYCR %>% 
  filter(Ks > 1.242339) %>% 
  filter(Ks < 3.159297) %>% 
  select(X)

write.table(C.crenata_CYCR_WGD_paralogs_peak1, quote = F, file = "Cyclopeltis_crenata_CYCR_WGD_paralogs_peak1.tsv")

#### Cyclopeltis presliana CYPR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.presliana_CYPR <- read.delim("ks_distributions/Cyclopeltis_presliana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.presliana_CYPR_filt <- C.presliana_CYPR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.presliana_CYPR_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.presliana_CYPR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cyclopeltis presliana (CYPR)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.presliana_CYPR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.presliana_CYPR_filt.0 <- C.presliana_CYPR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.presliana_CYPR_normalmixEM <- normalmixEM(C.presliana_CYPR_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peaks 
summary(C.presliana_CYPR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0968446 0.903155
#mu     0.1539342 2.177019
#sigma  0.1315557 0.954336
#loglik at estimate:  -8578.15  

# Generate list of paralogs 

C.presliana_CYPR_WGD_paralogs <- C.presliana_CYPR %>% 
  filter(Ks > 1.222683) %>% 
  filter(Ks < 3.131355) %>% 
  select(X)

write.table(C.presliana_CYPR_WGD_paralogs, quote = F, file = "Cyclopeltis_presliana_CYPR_WGD_paralogs.tsv")

#### Cyrtomium fortunei CYFO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.fortunei_CYFO <- read.delim("ks_distributions/Cyrtomium_fortunei.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.fortunei_CYFO_filt <- C.fortunei_CYFO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.fortunei_CYFO_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.fortunei_CYFO_filt, mapping = aes(x=Ks, ..scaled..*1300), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cyrtomium fortunei (CYFO)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.fortunei_CYFO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.fortunei_CYFO_filt.0 <- C.fortunei_CYFO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.fortunei_CYFO_normalmixEM <- normalmixEM(C.fortunei_CYFO_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peaks 
summary(C.fortunei_CYFO_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1134925 0.886507
#mu     0.1123298 2.114085
#sigma  0.0973522 0.978203
#loglik at estimate:  -8836.598   

# Generate list of paralogs 

C.fortunei_CYFO_WGD_paralogs <- C.fortunei_CYFO %>% 
  filter(Ks > 1.135882) %>% 
  filter(Ks < 3.092288) %>% 
  select(X)

write.table(C.fortunei_CYFO_WGD_paralogs, quote = F, file = "Cyrtomium_fortunei_CYFO_WGD_paralogs.tsv")

#### Cystopteris fragilis LHLE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.fragilis_LHLE <- read.delim("ks_distributions/Cystopteris_fragilis_LHLE.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.fragilis_LHLE_filt <- C.fragilis_LHLE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.fragilis_LHLE_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.fragilis_LHLE_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cystopteris fragilis (LHLE)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.fragilis_LHLE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.fragilis_LHLE_filt.0 <- C.fragilis_LHLE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.fragilis_LHLE_normalmixEM <- normalmixEM(C.fragilis_LHLE_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(C.fragilis_LHLE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0954135 0.904587
#mu     0.0503244 2.016453
#sigma  0.0437503 1.018637
#loglik at estimate:  -6394.741 

# Generate list of paralogs 

C.fragilis_LHLE_WGD_paralogs <- C.fragilis_LHLE %>% 
  filter(Ks > 0.997816) %>% 
  filter(Ks < 3.03509) %>% 
  select(X)

write.table(C.fragilis_LHLE_WGD_paralogs, quote = F, file = "Cystopteris_fragilis_LHLE_WGD_paralogs.tsv")

#### Cystopteris fragilis XXHP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.fragilis_XXHP <- read.delim("ks_distributions/Cystopteris_fragilis_XXHP.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.fragilis_XXHP_filt <- C.fragilis_XXHP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.fragilis_XXHP_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.fragilis_XXHP_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cystopteris fragilis_XXHP") + theme(plot.title = element_text(face = "italic"))

ggsave("C.fragilis_XXHP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.fragilis_XXHP_filt.0 <- C.fragilis_XXHP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.fragilis_XXHP_normalmixEM <- normalmixEM(C.fragilis_XXHP_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(C.fragilis_XXHP_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.131541 0.868459
#mu     0.124636 2.137332
#sigma  0.105723 0.991668
#loglik at estimate:  -7712.636  

# Generate list of paralogs 

C.fragilis_XXHP_WGD_paralogs <- C.fragilis_XXHP %>% 
  filter(Ks > 1.145664) %>% 
  filter(Ks < 3.129) %>% 
  select(X)

write.table(C.fragilis_XXHP_WGD_paralogs, quote = F, file = "Cystopteris_fragilis_XXHP_WGD_paralogs.tsv")

#### Cystopteris fragilis CYFR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.fragilis_CYFR <- read.delim("ks_distributions/Cystopteris_fragilis_CYFR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.fragilis_CYFR_filt <- C.fragilis_CYFR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.fragilis_CYFR_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.fragilis_CYFR_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cystopteris fragilis (CYFR)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.fragilis_CYFR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.fragilis_CYFR_filt.0 <- C.fragilis_CYFR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.fragilis_CYFR_normalmixEM <- normalmixEM(C.fragilis_CYFR_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(C.fragilis_CYFR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0727486 0.927251
#mu     0.1132628 2.061975
#sigma  0.0945926 0.987619
#loglik at estimate:  -6519.795 

# Generate list of paralogs 

C.fragilis_CYFR_WGD_paralogs <- C.fragilis_CYFR %>% 
  filter(Ks > 1.074356) %>% 
  filter(Ks < 3.049594) %>% 
  select(X)

write.table(C.fragilis_CYFR_WGD_paralogs, quote = F, file = "Cystopteris_fragilis_CYFR_WGD_paralogs.tsv")

#### Cystopteris reevesiana RICC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.reevesiana_RICC <- read.delim("ks_distributions/Cystopteris_reevesiana_RICC.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.reevesiana_RICC_filt <- C.reevesiana_RICC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.reevesiana_RICC_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.reevesiana_RICC_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cystopteris reevesiana (RICC)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.reevesiana_RICC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.reevesiana_RICC_filt.0 <- C.reevesiana_RICC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.reevesiana_RICC_normalmixEM <- normalmixEM(C.reevesiana_RICC_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(C.reevesiana_RICC_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.129480 0.870520
#mu     0.121651 2.079485
#sigma  0.101404 0.983465
#loglik at estimate:  -5129.681 

# Generate list of paralogs 

C.reevesiana_RICC_WGD_paralogs <- C.reevesiana_RICC %>% 
  filter(Ks > 1.09602) %>% 
  filter(Ks < 3.06295) %>% 
  select(X)

write.table(C.reevesiana_RICC_WGD_paralogs, quote = F, file = "Cystopteris_reevesiana_RICC_WGD_paralogs.tsv")

#### Cystopteris utahensis HNDZ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

C.utahensis_HNDZ <- read.delim("ks_distributions/Cystopteris_utahensis_HNDZ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

C.utahensis_HNDZ_filt <- C.utahensis_HNDZ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.reevesiana_RICC_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.reevesiana_RICC_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Cystopteris utahensis (HNDZ)") + theme(plot.title = element_text(face = "italic"))

ggsave("C.utahensis_HNDZ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

C.utahensis_HNDZ_filt.0 <- C.utahensis_HNDZ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

C.utahensis_HNDZ_normalmixEM <- normalmixEM(C.utahensis_HNDZ_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(C.utahensis_HNDZ_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1736935 0.826306
#mu     0.0351724 1.836412
#sigma  0.0288088 1.074646
#loglik at estimate:  -7848.864 

# Generate list of paralogs 

C.utahensis_HNDZ_WGD_paralogs <- C.utahensis_HNDZ %>% 
  filter(Ks > 0.761766) %>% 
  filter(Ks < 2.911058) %>% 
  select(X)

write.table(C.utahensis_HNDZ_WGD_paralogs, quote = F, file = "Cystopteris_utahensis_HNDZ_WGD_paralogs.tsv")

#### Danaea nodosa DANO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.nodosa_DANO <- read.delim("ks_distributions/Danaea_nodosa.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.nodosa_DANO_filt <- D.nodosa_DANO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.nodosa_DANO_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.nodosa_DANO_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Danaea nodosa (DANO)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.nodosa_DANO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.nodosa_DANO_filt.0 <- D.nodosa_DANO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.nodosa_DANO_normalmixEM <- normalmixEM(D.nodosa_DANO_filt.0$Ks, k=2)
# Fits to 2 WGD peaks 
summary(D.nodosa_DANO_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.630798 0.369202
#mu     1.122458 2.942077
#sigma  0.633342 0.623444
#loglik at estimate:  -8102.47 

# Generate list of paralogs 

D.nodosa_DANO_WGD_paralogs_peak1 <- D.nodosa_DANO %>% 
  filter(Ks > 0.489116) %>% 
  filter(Ks < 1.7558) %>% 
  select(X)

write.table(D.nodosa_DANO_WGD_paralogs_peak1, quote = F, file = "Danaea_nodosa_DANO_WGD_paralogs_peak1.tsv")

D.nodosa_DANO_WGD_paralogs_peak2 <- D.nodosa_DANO %>% 
  filter(Ks > 2.318633) %>% 
  filter(Ks < 3.565521) %>% 
  select(X)

write.table(D.nodosa_DANO_WGD_paralogs_peak2, quote = F, file = "Danaea_nodosa_DANO_WGD_paralogs_peak2.tsv")

#### Davallia griffithiana DAGR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.griffithiana_DAGR <- read.delim("ks_distributions/Davallia_griffithiana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.griffithiana_DAGR_filt <- D.griffithiana_DAGR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.griffithiana_DAGR_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.griffithiana_DAGR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Davallia griffithiana (DAGR)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.griffithiana_DAGR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.griffithiana_DAGR_filt.0 <- D.griffithiana_DAGR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.griffithiana_DAGR_normalmixEM <- normalmixEM(D.griffithiana_DAGR_filt.0$Ks, k=2)
# Fits to 2 WGD peaks 
summary(D.griffithiana_DAGR_normalmixEM )

#summary of normalmixEM object:
#comp 1   comp 2
#lambda 0.135886 0.864114
#mu     0.246056 2.247278
#sigma  0.189883 0.948714
#loglik at estimate:  -8411.715 

# Generate list of paralogs 

D.griffithiana_DAGR_WGD_paralogs <- D.griffithiana_DAGR %>% 
  filter(Ks > 1.298564) %>% 
  filter(Ks < 3.195992) %>% 
  select(X)

write.table(D.griffithiana_DAGR_WGD_paralogs, quote = F, file = "Davallia_griffithiana_DAGR_WGD_paralogs.tsv")

#### Davallia fejeensis OQWW ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.fejeensis_OQWW <- read.delim("ks_distributions/Davallia_fejeensis_QQWW.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.fejeensis_OQWW_filt <- D.fejeensis_OQWW %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.fejeensis_OQWW_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.fejeensis_OQWW_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Davallia fejeensis (QQWW)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.fejeensis_OQWW.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.fejeensis_OQWW_filt.0 <- D.fejeensis_OQWW_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

# No evidence of WGD 

#### Davallia repens DARE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.repens_DARE  <- read.delim("ks_distributions/Davallia_repens.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.repens_DARE_filt <- D.repens_DARE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(C.crenata_CYCR_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(C.crenata_CYCR_filt, mapping = aes(x=Ks, ..scaled..*400), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Davallia repens (DARE)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.repens_DARE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.repens_DARE_filt.0 <- D.repens_DARE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.repens_DARE_normalmixEM <- normalmixEM(D.repens_DARE_filt.0$Ks, k=2)
# Fits to 2 WGD peaks, comp 3 is overfit 
summary(D.repens_DARE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.183724 0.816276
#mu     0.283883 2.293131
#sigma  0.208215 0.957268
#loglik at estimate:  -5902.998 

# Generate list of paralogs 

D.repens_DARE_WGD_paralogs <- D.repens_DARE %>% 
  filter(Ks > 1.335863) %>% 
  filter(Ks < 3.250399) %>% 
  select(X)

write.table(D.repens_DARE_WGD_paralogs, quote = F, file = "Davallia_repens_DARE_WGD_paralogs.tsv")

#### Deparia lobato-crenata FCHS####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.lobato_crenata_FCHS <- read.delim("ks_distributions/Deparia_lobato-crenata_FCHS.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.lobato_crenata_FCHS_filt <- D.lobato_crenata_FCHS %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.lobato_crenata_FCHS_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.lobato_crenata_FCHS_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Deparia lobato-crenata (FCHS)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.lobato_crenata_FCHS.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.lobato_crenata_FCHS_filt.0 <- D.lobato_crenata_FCHS_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.lobato_crenata_FCHS_normalmixEM <- normalmixEM(D.lobato_crenata_FCHS_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(D.lobato_crenata_FCHS_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1291637 0.870836
#mu     0.1030746 1.974940
#sigma  0.0864552 0.982984
#loglik at estimate:  -4382.343 

# Generate list of paralogs 

D.lobato_crenata_FCHS_WGD_paralogs <- D.lobato_crenata_FCHS %>% 
  filter(Ks > 0.991956) %>% 
  filter(Ks < 2.957924) %>% 
  select(X)

write.table(D.lobato_crenata_FCHS_WGD_paralogs, quote = F, file = "Deparia_lobato_crenata_FCHS_WGD_paralogs.tsv")

#### Dennstaedtia davallioides MTGC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.davallioides_MTGC <- read.delim("ks_distributions/Dennstaedtia_davallioides_MTGC.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.davallioides_MTGC_filt <- D.davallioides_MTGC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.davallioides_MTGC_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.davallioides_MTGC_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dennstaedtia davallioides (MTGC)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.davallioides_MTGC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.davallioides_MTGC_filt.0 <- D.davallioides_MTGC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.davallioides_MTGC_normalmixEM <- normalmixEM(D.davallioides_MTGC_filt.0$Ks, k=2)
# Fits to 1 WGD peak  
summary(D.davallioides_MTGC_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0600443 0.939956
#mu     0.0521134 2.064547
#sigma  0.0435725 0.996947
#loglik at estimate:  -4318.32 

# Generate list of paralogs 

D.davallioides_MTGC_WGD_paralogs <- D.davallioides_MTGC %>% 
  filter(Ks > 0.991956) %>% 
  filter(Ks < 2.957924) %>% 
  select(X)

write.table(D.davallioides_MTGC_WGD_paralogs, quote = F, file = "Dennstaedtia_davallioides_MTGC_WGD_paralogs.tsv")

#### Dennstaedtia hirsuta DEHR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.hirsuta_DEHR <- read.delim("ks_distributions/Dennstaedtia_hirsuta_DEHR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.hirsuta_DEHR_filt <- D.hirsuta_DEHR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.hirsuta_DEHR_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.hirsuta_DEHR_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dennstaedtia hirsuta (DEHR)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.hirsuta_DEHR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.hirsuta_DEHR_filt.0 <- D.hirsuta_DEHR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.hirsuta_DEHR_normalmixEM <- normalmixEM(D.hirsuta_DEHR_filt.0$Ks, k=3)
# Fits to 1 WGD peak, comps 1 and 3 are overfit
summary(D.hirsuta_DEHR_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.0978466 0.673653 0.228501
#mu     0.1741911 1.800972 3.333123
#sigma  0.1429459 0.770088 0.400835
#loglik at estimate:  -9018.767 

# Generate list of paralogs 

D.hirsuta_DEHR_WGD_paralogs <- D.hirsuta_DEHR %>% 
  filter(Ks > 1.030884) %>% 
  filter(Ks < 2.57106) %>% 
  select(X)

write.table(D.hirsuta_DEHR_WGD_paralogs, quote = F, file = "Dennstaedtia_hirsuta_DEHR_WGD_paralogs.tsv")

#### Dennstaedtia hirsuta DEHI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.hirsuta_DEHI <- read.delim("ks_distributions/Dennstaedtia_hirsuta_DEHI.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.hirsuta_DEHI_filt <- D.hirsuta_DEHI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.hirsuta_DEHI_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.hirsuta_DEHI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dennstaedtia hirsuta (DEHI)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.hirsuta_DEHI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.hirsuta_DEHI_filt.0 <- D.hirsuta_DEHI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.hirsuta_DEHI_normalmixEM <- normalmixEM(D.hirsuta_DEHI_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comps 1 and 3 are overfit
summary(D.hirsuta_DEHI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0842787 0.915721
#mu     0.1197238 2.197849
#sigma  0.0964844 0.965090
#loglik at estimate:  -7824.506 

# Generate list of paralogs 

D.hirsuta_DEHI_WGD_paralogs <- D.hirsuta_DEHI %>% 
  filter(Ks > 1.232759) %>% 
  filter(Ks < 3.162939) %>% 
  select(X)

write.table(D.hirsuta_DEHI_WGD_paralogs, quote = F, file = "Dennstaedtia_hirsuta_DEHI_WGD_paralogs.tsv")

#### Dennstaedtia scabra DESC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.scabra_DESC <- read.delim("ks_distributions/Dennstaedtia_scabra.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.scabra_DESC_filt <- D.scabra_DESC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.scabra_DESC_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.scabra_DESC_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dennstaedtia scabra (DESC)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.scabra_DESC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.scabra_DESC_filt.0 <- D.scabra_DESC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.scabra_DESC_normalmixEM <- normalmixEM(D.scabra_DESC_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(D.scabra_DESC_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.101762 0.898238
#mu     0.144902 2.097377
#sigma  0.124479 0.983586
#loglik at estimate:  -8248.316

# Generate list of paralogs 

D.scabra_DESC_WGD_paralogs <- D.scabra_DESC %>% 
  filter(Ks > 1.113791) %>% 
  filter(Ks < 3.080963) %>% 
  select(X)

write.table(D.scabra_DESC_WGD_paralogs, quote = F, file = "Dennstaedtia_scabra_DESC_WGD_paralogs.tsv")

#### Deparia lancea DELA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.lancea_DELA <- read.delim("ks_distributions/Deparia_lancea.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.lancea_DELA_filt <- D.lancea_DELA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.lancea_DELA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.lancea_DELA_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Deparia lancea (DELA)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.lancea_DELA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.lancea_DELA_filt.0 <- D.lancea_DELA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.lancea_DELA_normalmixEM <- normalmixEM(D.lancea_DELA_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(D.lancea_DELA_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.866898 0.133102
#mu     2.031579 0.143847
#sigma  0.992617 0.107009
#loglik at estimate:  -9841.287 

# Generate list of paralogs 

D.lancea_DELA_WGD_paralogs <- D.lancea_DELA %>% 
  filter(Ks > 1.038962) %>% 
  filter(Ks < 3.024196) %>% 
  select(X)

write.table(D.lancea_DELA_WGD_paralogs, quote = F, file = "Deparia_lancea_DELA_WGD_paralogs.tsv")

#### Deparia okuboana DEOK ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.okuboana_DEOK <- read.delim("ks_distributions/Deparia_okuboana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.okuboana_DEOK_filt <- D.okuboana_DEOK %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.okuboana_DEOK_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.okuboana_DEOK_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Deparia okuboana (DEOK)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.okuboana_DEOK.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.okuboana_DEOK_filt.0 <- D.okuboana_DEOK_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.okuboana_DEOK_normalmixEM <- normalmixEM(D.okuboana_DEOK_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(D.okuboana_DEOK_normalmixEM )

#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.103342 0.896658
#mu     0.143006 2.107246
#sigma  0.118401 0.984763
#loglik at estimate:  -7646.756 

# Generate list of paralogs 

D.okuboana_DEOK_WGD_paralogs <- D.okuboana_DEOK %>% 
  filter(Ks > 1.122483) %>% 
  filter(Ks < 3.092009) %>% 
  select(X)

write.table(D.okuboana_DEOK_WGD_paralogs, quote = F, file = "Deparia_okuboana_DEOK_WGD_paralogs.tsv")

#### Deparia petersenii DEPE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.petersenii_DEPE <- read.delim("ks_distributions/Deparia_petersenii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.petersenii_DEPE_filt <- D.petersenii_DEPE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.petersenii_DEPE_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.petersenii_DEPE_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Deparia petersenii (DEPE)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.petersenii_DEPE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.petersenii_DEPE_filt.0 <- D.petersenii_DEPE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.petersenii_DEPE_normalmixEM <- normalmixEM(D.petersenii_DEPE_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(D.petersenii_DEPE_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.105628 0.894372
#mu     0.151134 2.042978
#sigma  0.102529 0.993361
#loglik at estimate:  -8085.352 

# Generate list of paralogs 

D.petersenii_DEPE_WGD_paralogs <- D.petersenii_DEPE %>% 
  filter(Ks > 1.049617) %>% 
  filter(Ks < 3.036339) %>% 
  select(X)

write.table(D.petersenii_DEPE_WGD_paralogs, quote = F, file = "Deparia_petersenii_DEPE_WGD_paralogs.tsv")

#### Dicksonia antarctica DIAT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.antarctica_DIAT <- read.delim("ks_distributions/Dicksonia_antarctica_DIAT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.antarctica_DIAT_filt <- D.antarctica_DIAT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.antarctica_DIAT_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.antarctica_DIAT_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dicksonia antarctica (DIAT)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.antarctica_DIAT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.antarctica_DIAT_filt.0 <- D.antarctica_DIAT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.antarctica_DIAT_normalmixEM <- normalmixEM(D.antarctica_DIAT_filt.0$Ks, k=5, maxit = 2000)
# Fits to 2 WGD peaks, comp 2 fits recent dups, comps 4 and 5 are overfit 
summary(D.antarctica_DIAT_normalmixEM )

#summary of normalmixEM object:
#          comp 1   comp 2   comp 3   comp 4    comp 5
#lambda 0.0162347 0.375340 0.259658 0.292621 0.0561462
#mu     0.5693677 0.242751 1.093225 2.366847 3.6156978
#sigma  0.0838104 0.121040 0.385228 0.628267 0.2249731
#loglik at estimate:  -8262.523 

# Generate list of paralogs 

D.antarctica_DIAT_WGD_paralogs_peak1 <- D.antarctica_DIAT %>% 
  filter(Ks > 0.4855573) %>% 
  filter(Ks < 0.6531781) %>% 
  select(X)

write.table(D.antarctica_DIAT_WGD_paralogs_peak1, quote = F, file = "Dicksonia_antarctica_DIAT_WGD_paralogs_peak1.tsv")

D.antarctica_DIAT_WGD_paralogs_peak2 <- D.antarctica_DIAT %>% 
  filter(Ks > 0.707997) %>% 
  filter(Ks < 1.478453) %>% 
  select(X)

write.table(D.antarctica_DIAT_WGD_paralogs_peak2, quote = F, file = "Dicksonia_antarctica_DIAT_WGD_paralogs_peak2.tsv")

#### Dicksonia antarctica DIAN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.antarctica_DIAN <- read.delim("ks_distributions/Dicksonia_antarctica_DIAN.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.antarctica_DIAN_filt <- D.antarctica_DIAN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.antarctica_DIAN_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.antarctica_DIAN_filt, mapping = aes(x=Ks, ..scaled..*1700), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dicksonia antarctica (DIAN)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.antarctica_DIAN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.antarctica_DIAN_filt.0 <- D.antarctica_DIAN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.antarctica_DIAN_normalmixEM <- normalmixEM(D.antarctica_DIAN_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, 3rd comp is overfit?
summary(D.antarctica_DIAN_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.374417 0.365919 0.259663
#mu     0.242382 1.226396 2.865266
#sigma  0.120624 0.530685 0.623430
#loglik at estimate:  -10107.79

# Generate list of paralogs 

D.antarctica_DIAN_WGD_paralogs_peak1 <- D.antarctica_DIAN %>% 
  filter(Ks > 0.121758) %>% 
  filter(Ks < 0.363006) %>% 
  select(X)

write.table(D.antarctica_DIAN_WGD_paralogs_peak1, quote = F, file = "Dicksonia_antarctica_DIAN_WGD_paralogs_peak1.tsv")

D.antarctica_DIAN_WGD_paralogs_peak2 <- D.antarctica_DIAN %>% 
  filter(Ks > 0.695711) %>% 
  filter(Ks < 1.757081) %>% 
  select(X)

write.table(D.antarctica_DIAN_WGD_paralogs_peak2, quote = F, file = "Dicksonia_antarctica_DIAN_WGD_paralogs_peak2.tsv")

#### Dicranopteris pedata DIPT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.pedata_DIPT <- read.delim("ks_distributions/Dicranopteris_pedata_DIPT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.pedata_DIPT_filt <- D.pedata_DIPT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.pedata_DIPT_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.pedata_DIPT_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dicranopteris pedata (DIPT)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.pedata_DIPT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.pedata_DIPT_filt.0 <- D.pedata_DIPT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.pedata_DIPT_normalmixEM <- normalmixEM(D.pedata_DIPT_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp2 is overfit? 
summary(D.pedata_DIPT_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.449388 0.550612
#mu     0.771653 2.509692
#sigma  0.428590 0.820831
#loglik at estimate:  -7114.401 

# Generate list of paralogs 

D.pedata_DIPT_WGD_paralogs <- D.pedata_DIPT %>% 
  filter(Ks > 0.343063) %>% 
  filter(Ks < 1.200243) %>% 
  select(X)

write.table(D.pedata_DIPT_WGD_paralogs, quote = F, file = "Dicranopteris_pedata_DIPT_WGD_paralogs_peak1.tsv")

#### Dicranopteris pedata DIPE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.pedata_DIPE <- read.delim("ks_distributions/Dicranopteris_pedata_DIPE.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.pedata_DIPE_filt <- D.pedata_DIPE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.pedata_DIPE_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.pedata_DIPE_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dicranopteris pedata (DIPE)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.pedata_DIPE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.pedata_DIPE_filt.0 <- D.pedata_DIPE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.pedata_DIPE_normalmixEM <- normalmixEM(D.pedata_DIPE_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp2 is overfit? 
summary(D.pedata_DIPE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.465248 0.534752
#mu     0.824363 2.550670
#sigma  0.427515 0.782986
#loglik at estimate:  -7987.194 

# Generate list of paralogs 

D.pedata_DIPE_WGD_paralogs <- D.pedata_DIPE %>% 
  filter(Ks > 0.396848) %>% 
  filter(Ks < 1.251878) %>% 
  select(X)

write.table(D.pedata_DIPE_WGD_paralogs, quote = F, file = "Dicranopteris_pedata_DIPE_WGD_paralogs_peak1.tsv")

#### Didymochlaena truncatula DITR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.truncatula_DITR <- read.delim("ks_distributions/Didymochlaena_truncatula.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.truncatula_DITR_filt <- D.truncatula_DITR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.truncatula_DITR_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.truncatula_DITR_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Didymochlaena truncatula (DITR)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.truncatula_DITR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.truncatula_DITR_filt.0 <- D.truncatula_DITR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.truncatula_DITR_normalmixEM <- normalmixEM(D.truncatula_DITR_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp2 is overfit? 
summary(D.truncatula_DITR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.108686 0.891314
#mu     0.128709 2.128453
#sigma  0.116182 0.977197
#loglik at estimate:  -7641.602

# Generate list of paralogs 

D.truncatula_DITR_WGD_paralogs <- D.truncatula_DITR %>% 
  filter(Ks > 1.151256) %>% 
  filter(Ks < 3.10565) %>% 
  select(X)

write.table(D.truncatula_DITR_WGD_paralogs, quote = F, file = "Didymochlaena_truncatula_DITR_WGD_paralogs_peak1.tsv")

#### Diplaziopsis cavaleriana DICA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.cavaleriana_DICA <- read.delim("ks_distributions/Diplaziopsis_cavaleriana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.cavaleriana_DICA_filt <- D.cavaleriana_DICA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.cavaleriana_DICA_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.cavaleriana_DICA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplaziopsis cavaleriana (DICA)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.cavaleriana_DICA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.cavaleriana_DICA_filt.0 <- D.cavaleriana_DICA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.cavaleriana_DICA_normalmixEM <- normalmixEM(D.cavaleriana_DICA_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp2 is overfit? 
summary(D.cavaleriana_DICA_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.770282 0.229718
#mu     1.528347 3.251008
#sigma  0.910012 0.446065
#loglik at estimate:  -8225.576 

# Generate list of paralogs 

D.cavaleriana_DICA_WGD_paralogs <- D.cavaleriana_DICA %>% 
  filter(Ks > 0.618335) %>% 
  filter(Ks < 2.438359) %>% 
  select(X)

write.table(D.cavaleriana_DICA_WGD_paralogs, quote = F, file = "Diplaziopsis_cavaleriana_DICA_WGD_paralogs.tsv")

#### Diplaziopsis javanica DIJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.javanica_DIJA <- read.delim("ks_distributions/Diplaziopsis_javanica.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.javanica_DIJA_filt <- D.javanica_DIJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.javanica_DIJA_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.javanica_DIJA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplaziopsis javanica (DIJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.javanica_DIJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.javanica_DIJA_filt.0 <- D.javanica_DIJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.javanica_DIJA_normalmixEM <- normalmixEM(D.javanica_DIJA_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak 
summary(D.javanica_DIJA_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1541005 0.845899
#mu     0.1234825 2.026886
#sigma  0.0952488 0.983156
#loglik at estimate:  -9615.373 

# Generate list of paralogs 

D.javanica_DIJA_WGD_paralogs <- D.javanica_DIJA %>% 
  filter(Ks > 1.04373) %>% 
  filter(Ks < 3.010042) %>% 
  select(X)

write.table(D.javanica_DIJA_WGD_paralogs, quote = F, file = "Diplaziopsis_javanica_DIJA_WGD_paralogs.tsv")

#### Diplazium chinense DICH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.chinense_DICH <- read.delim("ks_distributions/Diplazium_chinense.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.chinense_DICH_filt <- D.chinense_DICH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.chinense_DICH_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.chinense_DICH_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplazium chinense (DICH)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.chinense_DICH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.chinense_DICH_filt.0 <- D.chinense_DICH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.chinense_DICH_normalmixEM <- normalmixEM(D.chinense_DICH_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak 
summary(D.chinense_DICH_normalmixEM )
#summary of normalmixEM object:
#comp 1   comp 2
#lambda 9.49175e-05 0.999905
#mu     1.87257e+00 1.902067
#sigma  1.08394e+00 1.084478
#loglik at estimate:  -7158.177 

# Generate list of paralogs 

D.chinense_DICH_WGD_paralogs <- D.chinense_DICH %>% 
  filter(Ks > 0.817589) %>% 
  filter(Ks < 2.986545) %>% 
  select(X)

write.table(D.chinense_DICH_WGD_paralogs, quote = F, file = "Diplazium_chinense_DICH_WGD_paralogs.tsv")

#### Diplazium esculentum DIES ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.esculentum_DIES <- read.delim("ks_distributions/Diplazium_esculentum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.esculentum_DIES_filt <- D.esculentum_DIES %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.esculentum_DIES_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.esculentum_DIES_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplazium esculentum (DIES)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.esculentum_DIES.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD? 

#### Diplazium viridescens DIVI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.viridescens_DIVI <- read.delim("ks_distributions/Diplazium_viridescens.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.viridescens_DIVI_filt <- D.viridescens_DIVI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.viridescens_DIVI_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.viridescens_DIVI_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplazium viridescens (DIVI)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.viridescens_DIVI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.viridescens_DIVI_filt.0 <- D.viridescens_DIVI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.viridescens_DIVI_normalmixEM <- normalmixEM(D.viridescens_DIVI_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(D.viridescens_DIVI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.132846 0.867154
#mu     0.168704 2.055330
#sigma  0.132286 0.994305
#loglik at estimate:  -7886.918 

# Generate list of paralogs 

D.viridescens_DIVI_WGD_paralogs <- D.viridescens_DIVI %>% 
  filter(Ks > 1.061025) %>% 
  filter(Ks < 3.049635) %>% 
  select(X)

write.table(D.viridescens_DIVI_WGD_paralogs, quote = F, file = "Diplazium_viridescens_DIVI_WGD_paralogs.tsv")

#### Diplazium wichurae UFJN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.wichurae_UFJN <- read.delim("ks_distributions/Diplazium_wichurae_UFJN.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.wichurae_UFJN_filt <- D.wichurae_UFJN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.wichurae_UFJN_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.wichurae_UFJN_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplazium wichurae (UFJN)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.wichurae_UFJN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.wichurae_UFJN_filt.0 <- D.wichurae_UFJN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.wichurae_UFJN_normalmixEM <- normalmixEM(D.wichurae_UFJN_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(D.wichurae_UFJN_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1376308 0.862369
#mu     0.0987106 2.067924
#sigma  0.0848557 1.004542
#loglik at estimate:  -10211.91

# Generate list of paralogs 

D.wichurae_UFJN_WGD_paralogs <- D.wichurae_UFJN %>% 
  filter(Ks > 1.063382) %>% 
  filter(Ks < 3.072446) %>% 
  select(X)

write.table(D.wichurae_UFJN_WGD_paralogs, quote = F, file = "Diplazium_wichurae_UFJN_WGD_paralogs.tsv")

#### Diploterygium glucum DIGL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.glucum_DIGL <- read.delim("ks_distributions/Diplopterygium_glucum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.glucum_DIGL_filt <- D.glucum_DIGL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.glucum_DIGL_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.glucum_DIGL_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplopterygium glucum (DIGL)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.glucum_DIGL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.glucum_DIGL_filt.0 <- D.glucum_DIGL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.glucum_DIGL_normalmixEM <- normalmixEM(D.glucum_DIGL_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, 2nd comp is overfit? 
summary(D.glucum_DIGL_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2 
#lambda 0.480233 0.519767
#mu     0.648646 2.419961
#sigma  0.373309 0.837310
#loglik at estimate:  -8974.407

# Generate list of paralogs 

D.glucum_DIGL_WGD_paralogs <- D.glucum_DIGL %>% 
  filter(Ks > 0.275337) %>% 
  filter(Ks < 1.021955) %>% 
  select(X)

write.table(D.glucum_DIGL_WGD_paralogs, quote = F, file = "Diplopterygium_glucum_DIGL_WGD_paralogs.tsv")

#### Diploterygium laevissimum DILA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.laevissimum_DILA <- read.delim("ks_distributions/Diplopterygium_laevissimum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.laevissimum_DILA_filt <- D.laevissimum_DILA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.laevissimum_DILA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.laevissimum_DILA_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Diplopterygium laevissimum (DILA)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.laevissimum_DILA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.laevissimum_DILA_filt.0 <- D.laevissimum_DILA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.laevissimum_DILA_normalmixEM <- normalmixEM(D.laevissimum_DILA_filt.0$Ks, k=3, maxit = 2000)
# Fits to 2 WGD peaks,3rd comp is overfit 
summary(D.laevissimum_DILA_normalmixEM )

#summary of normalmixEM object:
#comp 1   comp 2   comp 3
#lambda 0.370096 0.511481 0.118423
#mu     0.674111 2.174627 3.511454
#sigma  0.358736 0.698376 0.289061
#loglik at estimate:  -7404.81 

# Generate list of paralogs 

D.laevissimum_DILA_WGD_paralogs_peak1 <- D.laevissimum_DILA %>% 
  filter(Ks > 0.315375) %>% 
  filter(Ks < 1.032847) %>% 
  select(X)

write.table(D.laevissimum_DILA_WGD_paralogs_peak1, quote = F, file = "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak1.tsv")

D.laevissimum_DILA_WGD_paralogs_peak2 <- D.laevissimum_DILA %>% 
  filter(Ks > 1.476251) %>% 
  filter(Ks < 2.873003) %>% 
  select(X)

write.table(D.laevissimum_DILA_WGD_paralogs_peak2, quote = F, file = "Diplopterygium_laevissimum_DILA_WGD_paralogs_peak2.tsv")

#### Dipteris conjugata DICO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.conjugata_DICO <- read.delim("ks_distributions/Dipteris_conjugata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.conjugata_DICO_filt <- D.conjugata_DICO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.conjugata_DICO_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.conjugata_DICO_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dipteris conjugata (DICO)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.conjugata_DICO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.conjugata_DICO_filt.0 <- D.conjugata_DICO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.conjugata_DICO_normalmixEM <- normalmixEM(D.conjugata_DICO_filt.0$Ks, k=2, maxit = 2000)
# Fits to 2 WGD peaks 
summary(D.conjugata_DICO_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.394323 0.605677
#mu     0.709739 2.516808
#sigma  0.427821 0.823381
#loglik at estimate:  -8161.053 

# Generate list of paralogs 

D.conjugata_DICO_WGD_paralogs_peak1 <- D.conjugata_DICO %>% 
  filter(Ks > 0.281918) %>% 
  filter(Ks < 1.13756) %>% 
  select(X)

write.table(D.conjugata_DICO_WGD_paralogs_peak1, quote = F, file = "Dipteris_conjugata_DICO_WGD_paralogs_peak1.tsv")

D.conjugata_DICO_WGD_paralogs_peak2 <- D.conjugata_DICO %>% 
  filter(Ks > 1.693427) %>% 
  filter(Ks < 3.340189) %>% 
  select(X)

write.table(D.conjugata_DICO_WGD_paralogs_peak2, quote = F, file = "Dipteris_conjugata_DICO_WGD_paralogs_peak2.tsv")

#### Dipteris lobiana DILO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.lobiana_DILO <- read.delim("ks_distributions/Dipteris_lobbiana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.lobiana_DILO_filt <- D.lobiana_DILO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.lobiana_DILO_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.lobiana_DILO_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dipteris lobiana (DILO)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.lobiana_DILO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.lobiana_DILO_filt.0 <- D.lobiana_DILO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.lobiana_DILO_normalmixEM <- normalmixEM(D.lobiana_DILO_filt.0$Ks, k=3, maxit = 2000)
# Fits to 2 WGD peaks, comp 3 is overfit  
summary(D.lobiana_DILO_normalmixEM )

#summary of normalmixEM object:
#         comp 1    comp 2   comp 3
#lambda 0.249770 0.2971246 0.453105
#mu     0.739047 0.1482494 2.407783
#sigma  0.334672 0.0934289 0.875619
#loglik at estimate:  -8845.068 

# Generate list of paralogs 

D.lobiana_DILO_WGD_paralogs_peak1 <- D.lobiana_DILO %>% 
  filter(Ks > 0.0548205) %>% 
  filter(Ks < 0.2416783) %>% 
  select(X)

write.table(D.lobiana_DILO_WGD_paralogs_peak1, quote = F, file = "Dipteris_lobiana_DILO_WGD_paralogs_peak1.tsv")

D.lobiana_DILO_WGD_paralogs_peak2 <- D.lobiana_DILO %>% 
  filter(Ks > 0.404375) %>% 
  filter(Ks < 1.073719) %>% 
  select(X)

write.table(D.lobiana_DILO_WGD_paralogs_peak2, quote = F, file = "Dipteris_lobiana_DILO_WGD_paralogs_peak2.tsv")

#### Dryopteris decipiens DRDE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.decipiens_DRDE <- read.delim("ks_distributions/Dryopteris_decipiens.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.decipiens_DRDE_filt <- D.decipiens_DRDE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.decipiens_DRDE_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.decipiens_DRDE_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dryopteris decipiens (DRDE)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.decipiens_DRDE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.decipiens_DRDE_filt.0 <- D.decipiens_DRDE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.decipiens_DRDE_normalmixEM <- normalmixEM(D.decipiens_DRDE_filt.0$Ks, k=2, maxit = 2000)
# Fits to 2 WGD peaks, comp 3 is overfit  
summary(D.decipiens_DRDE_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.164435 0.835565
#mu     0.151208 2.047509
#sigma  0.111468 0.992116
#loglik at estimate:  -9579.876

# Generate list of paralogs 

D.decipiens_DRDE_WGD_paralogs <- D.decipiens_DRDE %>% 
  filter(Ks > 1.055393) %>% 
  filter(Ks < 3.039625) %>% 
  select(X)

write.table(D.decipiens_DRDE_WGD_paralogs, quote = F, file = "Dryopteris_decipiens_DRDE_WGD_paralogs.tsv")

#### Dryopteris pseudocaenopteris DRPS ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

D.pseudocaenopteris_DRPS <- read.delim("ks_distributions/Dryopteris_pseudocaenopteris.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

D.pseudocaenopteris_DRPS_filt <- D.pseudocaenopteris_DRPS %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(D.pseudocaenopteris_DRPS_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(D.pseudocaenopteris_DRPS_filt, mapping = aes(x=Ks, ..scaled..*400), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Dryopteris pseudocaenopteris (DRPS)") + theme(plot.title = element_text(face = "italic"))

ggsave("D.pseudocaenopteris_DRPS.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

D.pseudocaenopteris_DRPS_filt.0 <- D.pseudocaenopteris_DRPS_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

D.pseudocaenopteris_DRPS_normalmixEM <- normalmixEM(D.pseudocaenopteris_DRPS_filt.0$Ks, k=3, maxit = 2000)
# Fits to 1 WGD peak, comp 3 is overfit
summary(D.pseudocaenopteris_DRPS_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.103058 0.666394 0.230548
#mu     0.137803 1.734409 3.316894
#sigma  0.110004 0.757696 0.394747
#loglik at estimate:  -6795.518  

# Generate list of paralogs 

D.pseudocaenopteris_DRPS_WGD_paralogs <- D.pseudocaenopteris_DRPS %>% 
  filter(Ks > 0.976713) %>% 
  filter(Ks < 2.492105) %>% 
  select(X)

write.table(D.pseudocaenopteris_DRPS_WGD_paralogs, quote = F, file = "Dryopteris_pseudocaenopteris_DRPS_WGD_paralogs.tsv")

#### Elaphoglossum mcclurei ELMC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.mcclurei_ELMC <- read.delim("ks_distributions/Elaphoglossum_mccluri.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.mcclurei_ELMC_filt <- E.mcclurei_ELMC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.mcclurei_ELMC_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.mcclurei_ELMC_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Elaphoglossum mcclurei (ELMC)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.mcclurei_ELMC.png", height = 5, width = 8, dpi = 300)

# No evidence for WGD 

#### Elaphoglossum yoshinagae ELYO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.yoshinagae_ELYO <- read.delim("ks_distributions/Elaphoglossum_yoshinagae.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.yoshinagae_ELYO_filt <- E.yoshinagae_ELYO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.yoshinagae_ELYO_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.yoshinagae_ELYO_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Elaphoglossum yoshinagae (ELYO)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.yoshinagae_ELYO.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 
 
#### Equisetum arvense EQAR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.arvense_EQAR <- read.delim("ks_distributions/Equisetum_arvense.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.arvense_EQAR_filt <- E.arvense_EQAR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.arvense_EQAR_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.arvense_EQAR_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Equisetum arvense (EQAR)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.arvense_EQAR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

E.arvense_EQAR_filt.0 <- E.arvense_EQAR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

E.arvense_EQAR_normalmixEM <- normalmixEM(E.arvense_EQAR_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(E.arvense_EQAR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.587800 0.412200
#mu     2.411727 0.783972
#sigma  0.834738 0.446547
#loglik at estimate:  -8162.632 

# Generate list of paralogs 

E.arvense_EQAR_WGD_paralogs <- E.arvense_EQAR %>% 
  filter(Ks > 0.337425) %>% 
  filter(Ks < 1.230519) %>% 
  select(X)

write.table(E.arvense_EQAR_WGD_paralogs, quote = F, file = "Equisetum_arvense_EQAR_WGD_paralogs.tsv")

#### Equisetum diffusum EQDI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.diffusum_EQDI <- read.delim("ks_distributions/Equisetum_diffusum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.diffusum_EQDI_filt <- E.diffusum_EQDI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.diffusum_EQDI_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.diffusum_EQDI_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Equisetum diffusum (EQDI)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.diffusum_EQDI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

E.diffusum_EQDI_filt.0 <- E.diffusum_EQDI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

E.diffusum_EQDI_normalmixEM <- normalmixEM(E.diffusum_EQDI_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(E.diffusum_EQDI_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.469935 0.530065
#mu     0.729876 2.366827
#sigma  0.438441 0.869361
#loglik at estimate:  -8608.775

# Generate list of paralogs 

E.diffusum_EQDI_WGD_paralogs <- E.diffusum_EQDI %>% 
  filter(Ks > 0.291435) %>% 
  filter(Ks < 1.168317) %>% 
  select(X)

write.table(E.diffusum_EQDI_WGD_paralogs, quote = F, file = "Equisetum_diffusum_EQDI_WGD_paralogs.tsv")

#### Equisetum diffusum CAPN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.diffusum_CAPN <- read.delim("ks_distributions/Equisetum_diffusum_CAPN.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.diffusum_CAPN_filt <- E.diffusum_CAPN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.diffusum_CAPN_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.diffusum_CAPN_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Equisetum diffusum (CAPN)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.diffusum_CAPN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

E.diffusum_CAPN_filt.0 <- E.diffusum_CAPN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

E.diffusum_CAPN_normalmixEM <- normalmixEM(E.diffusum_CAPN_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(E.diffusum_CAPN_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.407691 0.592309
#mu     0.674625 2.345005
#sigma  0.439870 0.851709
#loglik at estimate:  -6107.819 

# Generate list of paralogs 

E.diffusum_CAPN_WGD_paralogs <- E.diffusum_CAPN %>% 
  filter(Ks > 0.234755) %>% 
  filter(Ks < 1.114495) %>% 
  select(X)

write.table(E.diffusum_CAPN_WGD_paralogs, quote = F, file = "Equisetum_diffusum_CAPN_WGD_paralogs.tsv")

#### Equisetum hyemale JVSZ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

E.hymale_JVSZ <- read.delim("ks_distributions/Equisetum_hyemale_JVSZ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

E.hymale_JVSZ_filt <- E.hymale_JVSZ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(E.hymale_JVSZ_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(E.hymale_JVSZ_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Equisetum hyemale (JVSZ)") + theme(plot.title = element_text(face = "italic"))

ggsave("E.hymale_JVSZ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

E.hymale_JVSZ_filt.0 <- E.hymale_JVSZ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

E.hymale_JVSZ_normalmixEM <- normalmixEM(E.hymale_JVSZ_filt.0$Ks, k=2)
# Fits to 1 WGD peak- first component is the peak that is likely real 
summary(E.hymale_JVSZ_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.416046 0.583954
#mu     0.736772 2.327417
#sigma  0.439732 0.873584
#loglik at estimate:  -10656.83 

# Generate list of paralogs 

E.hymale_JVSZ_WGD_paralogs <- E.hymale_JVSZ %>% 
  filter(Ks > 0.29704) %>% 
  filter(Ks < 1.176504) %>% 
  select(X)

write.table(E.hymale_JVSZ_WGD_paralogs, quote = F, file = "Equisetum_hymale_JVSZ_WGD_paralogs.tsv")

#### Gaga arizonica DCDT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.arizonica_DCDT <- read.delim("ks_distributions/Gaga_arizonica_DCDT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.arizonica_DCDT_filt <- G.arizonica_DCDT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.arizonica_DCDT_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.arizonica_DCDT_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Gaga arizonica (DCDT)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.arizonica_DCDT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

G.arizonica_DCDT_filt.0 <- G.arizonica_DCDT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

G.arizonica_DCDT_normalmixEM <- normalmixEM(G.arizonica_DCDT_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(G.arizonica_DCDT_normalmixEM )

#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.141301 0.858699
#mu     0.156991 2.271419
#sigma  0.139954 0.947753
#loglik at estimate:  -5733.549 

# Generate list of paralogs 

G.arizonica_DCDT_WGD_paralogs <- G.arizonica_DCDT %>% 
  filter(Ks > 1.323666) %>% 
  filter(Ks < 3.219172) %>% 
  select(X)

write.table(G.arizonica_DCDT_WGD_paralogs, quote = F, file = "Gaga_arizonica_DCDT_WGD_paralogs.tsv")

#### Goniophlebium niponicum GONP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.niponicum_GONP <- read.delim("ks_distributions/Goniophlebium_niponicum_GONP.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.niponicum_GONP_filt <- G.niponicum_GONP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.niponicum_GONP_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.niponicum_GONP_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Goniophlebium niponicum (GONP)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.niponicum_GONP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

G.niponicum_GONP_filt.0 <- G.niponicum_GONP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

G.niponicum_GONP_normalmixEM <- normalmixEM(G.niponicum_GONP_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(G.niponicum_GONP_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.147583 0.852417
#mu     0.244697 2.313420
#sigma  0.187435 0.929015
#loglik at estimate:  -8535.079

# Generate list of paralogs 

G.niponicum_GONP_WGD_paralogs <- G.niponicum_GONP %>% 
  filter(Ks > 1.384405) %>% 
  filter(Ks < 3.242435) %>% 
  select(X)

write.table(G.niponicum_GONP_WGD_paralogs, quote = F, file = "Goniophlebium_niponicum_GONP_WGD_paralogs.tsv")

#### Goniophlebium niponicum GONI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.niponicum_GONI <- read.delim("ks_distributions/Goniophlebium_niponicum_GONI.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.niponicum_GONI_filt <- G.niponicum_GONI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.niponicum_GONI_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.niponicum_GONI_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Goniophlebium niponicum (GONI)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.niponicum_GONI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

G.niponicum_GONI_filt.0 <- G.niponicum_GONI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

G.niponicum_GONI_normalmixEM <- normalmixEM(G.niponicum_GONI_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(G.niponicum_GONI_normalmixEM )

#summary of normalmixEM object:
#comp 1   comp 2
#lambda 0.127135 0.872865
#mu     0.177297 2.213161
#sigma  0.118665 1.018220
#loglik at estimate:  -8544.997 

# Generate list of paralogs 

G.niponicum_GONI_WGD_paralogs <- G.niponicum_GONI %>% 
  filter(Ks > 1.194941) %>% 
  filter(Ks < 3.231381) %>% 
  select(X)

write.table(G.niponicum_GONI_WGD_paralogs, quote = F, file = "Goniophlebium_niponicum_GONI_WGD_paralogs.tsv")

#### Grammitis dorsipila GRDO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.dorsipila_GORD <- read.delim("ks_distributions/Grammitis_dorsipila.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.dorsipila_GORD_filt <- G.dorsipila_GORD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.dorsipila_GORD_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.dorsipila_GORD_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Grammitis dorsipila (GRDO)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.dorsipila_GRDO.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Gymnocarpium dryopteris HEGQ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.dryopteris_HEGQ <- read.delim("ks_distributions/Gymnocarpium_dryopteris_HEGQ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.dryopteris_HEGQ_filt <- G.dryopteris_HEGQ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.dryopteris_HEGQ_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.dryopteris_HEGQ_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Gymnocarpium dryopteris (HEGQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.dryopteris_HEGQ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

G.dryopteris_HEGQ_filt.0 <- G.dryopteris_HEGQ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

G.dryopteris_HEGQ_normalmixEM <- normalmixEM(G.dryopteris_HEGQ_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(G.dryopteris_HEGQ_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1386472 0.861353
#mu     0.0949137 2.020470
#sigma  0.0851333 1.017617
#loglik at estimate:  -6525.74 

# Generate list of paralogs 

G.dryopteris_HEGQ_WGD_paralogs <- G.dryopteris_HEGQ %>% 
  filter(Ks > 1.002853) %>% 
  filter(Ks < 3.038087) %>% 
  select(X)

write.table(G.dryopteris_HEGQ_WGD_paralogs, quote = F, file = "Gymnocarpium_dryopteris_HEGQ_WGD_paralogs.tsv")

#### Gymnocarpium oyamense GYOY ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

G.oyamense_GYOY <- read.delim("ks_distributions/Gymnocarpium_oyamense.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

G.oyamense_GYOY_filt <- G.oyamense_GYOY %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(G.oyamense_GYOY_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(G.oyamense_GYOY_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Gymnocarpium oyamense (GYOY)") + theme(plot.title = element_text(face = "italic"))

ggsave("G.oyamense_GYOY.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

G.oyamense_GYOY_filt.0 <- G.oyamense_GYOY_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

G.oyamense_GYOY_normalmixEM <- normalmixEM(G.oyamense_GYOY_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(G.oyamense_GYOY_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.116811 0.883189
#mu     0.194597 2.077469
#sigma  0.157470 0.988081
#loglik at estimate:  -9556.057 

# Generate list of paralogs 

G.oyamense_GYOY_WGD_paralogs <-G.oyamense_GYOY %>% 
  filter(Ks > 1.089388) %>% 
  filter(Ks < 3.06555) %>% 
  select(X)

write.table(G.oyamense_GYOY_WGD_paralogs, quote = F, file = "Gymnocarpium_oyamense_GYOY_WGD_paralogs.tsv")

#### Haplopteris amboinensis HAAM ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.amboinensis_HAAM <- read.delim("ks_distributions/Haplopteris_amboinensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.amboinensis_HAAM_filt <- H.amboinensis_HAAM %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.amboinensis_HAAM_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.amboinensis_HAAM_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Haplopteris amboinensis (HAAM)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.amboinensis_HAAM.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.amboinensis_HAAM_filt.0 <- H.amboinensis_HAAM_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.amboinensis_HAAM_normalmixEM <- normalmixEM(H.amboinensis_HAAM_filt.0$Ks, k=2)
# Fits to 2 WGD peak
summary(H.amboinensis_HAAM_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.420165 0.579835
#mu     0.475503 2.402147
#sigma  0.240105 0.898828
#loglik at estimate:  -6861.701 

# Generate list of paralogs 

H.amboinensis_HAAM_WGD_paralogs_peak1 <-H.amboinensis_HAAM %>% 
  filter(Ks > 0.235398) %>% 
  filter(Ks < 0.715608) %>% 
  select(X)

write.table(H.amboinensis_HAAM_WGD_paralogs_peak1, quote = F, file = "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak1.tsv")

H.amboinensis_HAAM_WGD_paralogs_peak2 <-H.amboinensis_HAAM %>% 
  filter(Ks > 1.503319) %>% 
  filter(Ks < 3.300975) %>% 
  select(X)

write.table(H.amboinensis_HAAM_WGD_paralogs_peak2, quote = F, file = "Haplopteris_amboinensis_HAAM_WGD_paralogs_peak2.tsv")

#### Haplopteris heterophylla HAHE #####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.heterophylla_HAHE <- read.delim("ks_distributions/Haplopteris_heterophylla.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.heterophylla_HAHE_filt <- H.heterophylla_HAHE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.heterophylla_HAHE_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.heterophylla_HAHE_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Haplopteris heterophylla (HAHE)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.heterophylla_HAHE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.heterophylla_HAHE_filt.0 <- H.heterophylla_HAHE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.heterophylla_HAHE_normalmixEM <- normalmixEM(H.heterophylla_HAHE_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.heterophylla_HAHE_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.372718 0.627282
#mu     0.515454 2.409398
#sigma  0.288848 0.887838
#loglik at estimate:  -8115.113   

# Generate list of paralogs 

H.heterophylla_HAHE_WGD_paralogs_peak1 <-H.heterophylla_HAHE %>% 
  filter(Ks > 0.226606) %>% 
  filter(Ks < 0.804302) %>% 
  select(X)

write.table(H.heterophylla_HAHE_WGD_paralogs_peak1, quote = F, file = "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak1.tsv")

H.heterophylla_HAHE_WGD_paralogs_peak2 <-H.heterophylla_HAHE %>% 
  filter(Ks > 1.52156) %>% 
  filter(Ks < 3.297236) %>% 
  select(X)

write.table(H.heterophylla_HAHE_WGD_paralogs_peak2, quote = F, file = "Haplopteris_heterophylla_HAHE_WGD_paralogs_peak2.tsv")

#### Hemionitis arifolia HEAR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.arifolia_HEAR <- read.delim("ks_distributions/Heminoitis_arifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.arifolia_HEAR_filt <- H.arifolia_HEAR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.arifolia_HEAR_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.arifolia_HEAR_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Heminoitis arifolia (HEAR)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.arifolia_HEAR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.arifolia_HEAR_filt.0 <- H.arifolia_HEAR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.arifolia_HEAR_normalmixEM <- normalmixEM(H.arifolia_HEAR_filt.0$Ks, k=2)
# Fits to 2 WGD peak
summary(H.arifolia_HEAR_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1068890 0.893111
#mu     0.1043397 2.284217
#sigma  0.0940336 0.955953
#loglik at estimate:  -5655.187 

# Generate list of paralogs 

H.arifolia_HEAR_WGD_paralogs <-H.arifolia_HEAR %>% 
  filter(Ks > 1.328264) %>% 
  filter(Ks < 3.24017) %>% 
  select(X)

write.table(H.arifolia_HEAR_WGD_paralogs, quote = F, file = "Heminoitis_arifolia_HEAR_HAHE_WGD_paralogs.tsv")

#### Histiopteris incisa HIIC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.incisa_HIIC <- read.delim("ks_distributions/Histiopteris_incisa_HIIC.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.incisa_HIIC_filt <- H.incisa_HIIC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.incisa_HIIC_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.incisa_HIIC_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Histiopteris incisa (HIIC)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.incisa_HIIC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.incisa_HIIC_filt.0 <- H.incisa_HIIC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.incisa_HIIC_normalmixEM <- normalmixEM(H.incisa_HIIC_filt.0$Ks, k=3)
# Fits to 1 WGD, ignore 3rd comp
summary(H.incisa_HIIC_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.100891 0.743143 0.155966
#mu     0.166498 1.856791 3.411350
#sigma  0.128334 0.789016 0.348935
#loglik at estimate:  -7920.677

# Generate list of paralogs 

H.incisa_HIIC_WGD_paralogs <-H.incisa_HIIC %>% 
  filter(Ks > 1.067775) %>% 
  filter(Ks < 2.645807) %>% 
  select(X)

write.table(H.incisa_HIIC_WGD_paralogs, quote = F, file = "Histiopteris_incisa_HIIC_WGD_paralogs.tsv")

#### Histiopteris incisa HIIN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.incisa_HIIN <- read.delim("ks_distributions/Histiopteris_incisa_HIIN.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.incisa_HIIN_filt <- H.incisa_HIIN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.incisa_HIIN_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.incisa_HIIN_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Histiopteris incisa (HIIN)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.incisa_HIIN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.incisa_HIIN_filt.0 <- H.incisa_HIIN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.incisa_HIIN_normalmixEM <- normalmixEM(H.incisa_HIIN_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.incisa_HIIN_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.144096 0.855904
#mu     0.169212 2.087649
#sigma  0.132722 0.964446
#loglik at estimate:  -5451.221 

# Generate list of paralogs 

H.incisa_HIIN_WGD_paralogs <-H.incisa_HIIN %>% 
  filter(Ks > 1.123203) %>% 
  filter(Ks < 3.052095) %>% 
  select(X)

write.table(H.incisa_HIIN_WGD_paralogs, quote = F, file = "Histiopteris_incisa_HIIN_WGD_paralogs.tsv")

#### Homalosorus pycnocarpos OCZL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.pycnocarpos_OCZL <- read.delim("ks_distributions/Homalosorus_pycnocarpos_OCZL.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.pycnocarpos_OCZL_filt <- H.pycnocarpos_OCZL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.pycnocarpos_OCZL_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.pycnocarpos_OCZL_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Homalosorus pycnocarpos (OCZL)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.pycnocarpos_OCZL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.pycnocarpos_OCZL_filt.0 <- H.pycnocarpos_OCZL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.pycnocarpos_OCZL_normalmixEM <- normalmixEM(H.pycnocarpos_OCZL_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(H.pycnocarpos_OCZL_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.124394 0.875606
#mu     0.138112 2.088440
#sigma  0.103920 0.974306
#loglik at estimate:  -5821.859

# Generate list of paralogs 

H.pycnocarpos_OCZL_WGD_paralogs <-H.pycnocarpos_OCZL %>% 
  filter(Ks > 1.114134) %>% 
  filter(Ks < 3.062746) %>% 
  select(X)

write.table(H.pycnocarpos_OCZL_WGD_paralogs, quote = F, file = "Homalosorus_pycnocarpos_OCZL_WGD_paralogs.tsv")

#### Hymenasplenium sp. HYXQ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.sp_2018 <- read.delim("ks_distributions/Hymenasplenium_sp._XQ-2018.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.sp_2018_filt <- H.sp_2018 %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.sp_2018_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.sp_2018_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hymenasplenium sp. (HYXQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.sp_2018_HYXQ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.sp_2018_filt.0 <- H.sp_2018_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.sp_2018_normalmixEM <- normalmixEM(H.sp_2018_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.sp_2018_normalmixEM )

#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.140341 0.859659
#mu     0.149089 2.115849
#sigma  0.104989 1.015158
#loglik at estimate:  -8485.043

# Generate list of paralogs 

H.sp_2018_WGD_paralogs <-H.sp_2018 %>% 
  filter(Ks > 1.100691) %>% 
  filter(Ks < 3.131007) %>% 
  select(X)

write.table(H.sp_2018_WGD_paralogs, quote = F, file = "Hymenasplenium_sp_HYXQ_WGD_paralogs.tsv")

#### Hymenophyllum holochilum HYHO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.holochilum_HYHO <- read.delim("ks_distributions/Hymenophyllum_holochium.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.holochilum_HYHO_filt <- H.holochilum_HYHO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.holochilum_HYHO_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.holochilum_HYHO_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hymenophyllum holochilum (HYHO)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.holochilum_HYHO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.holochilum_HYHO_filt.0 <- H.holochilum_HYHO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.holochilum_HYHO_normalmixEM <- normalmixEM(H.holochilum_HYHO_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(H.holochilum_HYHO_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.412524 0.587476
#mu     0.669666 2.559384
#sigma  0.433643 0.815032
#loglik at estimate:  -7090.335 

# Generate list of paralogs 

H.holochilum_HYHO_WGD_paralogs <-H.holochilum_HYHO %>% 
  filter(Ks > 0.236023) %>% 
  filter(Ks < 1.103309) %>% 
  select(X)

write.table(H.holochilum_HYHO_WGD_paralogs, quote = F, file = "Hymenophyllum_holochilum_HYHO_WGD_paralogs.tsv")

#### Hymenophyllum sp. HYME ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.sp_HYME <- read.delim("ks_distributions/Hymenophyllum_sp._XQ-2018.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.sp_HYME_filt <- H.sp_HYME %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.sp_HYME_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.sp_HYME_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hymenophyllum sp. (HYME)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.sp_HYME.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.sp_HYME_filt.0 <- H.sp_HYME_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.sp_HYME_normalmixEM <- normalmixEM(H.sp_HYME_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.sp_HYME_normalmixEM )

#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.328772 0.671228
#mu     0.697396 2.357022
#sigma  0.458564 0.862145
#loglik at estimate:  -9666.23  

# Generate list of paralogs 

H.sp_HYME_WGD_paralogs <-H.sp_HYME %>% 
  filter(Ks > 1.494877) %>% 
  filter(Ks < 3.219167) %>% 
  select(X)

write.table(H.sp_HYME_WGD_paralogs, quote = F, file = "Hymenophyllum_sp_HYME_WGD_paralogs.tsv")

#### Hypodematium crenatum HYCT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.crenatum_HYCT <- read.delim("ks_distributions/Hypodematium_crenatum_HYCT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.crenatum_HYCT_filt <- H.crenatum_HYCT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.crenatum_HYCT_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.crenatum_HYCT_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hypodematium crenatum (HYCT)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.crenatum_HYCT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.crenatum_HYCT_filt.0 <- H.crenatum_HYCT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.crenatum_HYCT_normalmixEM <- normalmixEM(H.crenatum_HYCT_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(H.crenatum_HYCT_normalmixEM )
#summary of normalmixEM object:
#  comp 1    comp 2
#lambda 0.922187 0.0778126
#mu     2.201535 0.1589683
#sigma  0.954867 0.1203547
#loglik at estimate:  -7352.714 

# Generate list of paralogs 

H.crenatum_HYCT_WGD_paralogs <-H.crenatum_HYCT %>% 
  filter(Ks > 1.246668) %>% 
  filter(Ks < 3.156402) %>% 
  select(X)

write.table(H.crenatum_HYCT_WGD_paralogs, quote = F, file = "Hypodematium_crenatum_HYCT_WGD_paralogs.tsv")

#### Hypodematium crenatum HYCR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.crenatum_HYCR <- read.delim("ks_distributions/Hypodematium_crenatum_HYCR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.crenatum_HYCR_filt <- H.crenatum_HYCR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.crenatum_HYCR_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.crenatum_HYCR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hypodematium crenatum (HYCR)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.crenatum_HYCR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.crenatum_HYCR_filt.0 <- H.crenatum_HYCR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.crenatum_HYCR_normalmixEM <- normalmixEM(H.crenatum_HYCR_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.crenatum_HYCR_normalmixEM )

#summary of normalmixEM object:
#comp 1   comp 2
#lambda 0.0859767 0.914023
#mu     0.1010347 2.111886
#sigma  0.0933855 0.988843
#loglik at estimate:  -7855.809

# Generate list of paralogs 

H.crenatum_HYCR_WGD_paralogs <-H.crenatum_HYCR %>% 
  filter(Ks > 1.123043) %>% 
  filter(Ks < 3.100729) %>% 
  select(X)

write.table(H.crenatum_HYCR_WGD_paralogs, quote = F, file = "Hypodematium_crenatum_HYCR_WGD_paralogs.tsv")

#### Hypolepis punctata HYPT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.punctata_HYPT <- read.delim("ks_distributions/Hypolepis_punctata_HYPT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.punctata_HYPT_filt <- H.punctata_HYPT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.punctata_HYPT_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.punctata_HYPT_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hypolepis punctata (HYPT)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.punctata_HYPT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.punctata_HYPT_filt.0 <- H.punctata_HYPT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.punctata_HYPT_normalmixEM <- normalmixEM(H.punctata_HYPT_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(H.punctata_HYPT_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.104275 0.895725
#mu     0.116483 2.058803
#sigma  0.105046 1.011626
#loglik at estimate:  -7473.087 

# Generate list of paralogs 

H.punctata_HYPT_WGD_paralogs <-H.punctata_HYPT %>% 
  filter(Ks > 1.047177) %>% 
  filter(Ks < 3.070429) %>% 
  select(X)

write.table(H.punctata_HYPT_WGD_paralogs, quote = F, file = "Hypolepis_punctata_HYPT_WGD_paralogs.tsv")

#### Hypolepis punctata HYPU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

H.punctata_HYPU <- read.delim("ks_distributions/Hypolepis_punctata_HYPU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

H.punctata_HYPU_filt <- H.punctata_HYPU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(H.punctata_HYPU_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(H.punctata_HYPU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Hypolepis punctata (HYPU)") + theme(plot.title = element_text(face = "italic"))

ggsave("H.punctata_HYPU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

H.punctata_HYPU_filt.0 <- H.punctata_HYPU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

H.punctata_HYPU_normalmixEM <- normalmixEM(H.punctata_HYPU_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(H.punctata_HYPU_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.124381 0.875619
#mu     0.106011 2.019084
#sigma  0.089678 0.987376
#loglik at estimate:  -7377.683  

# Generate list of paralogs 

H.punctata_HYPU_WGD_paralogs <-H.punctata_HYPU %>% 
  filter(Ks > 1.031708) %>% 
  filter(Ks < 3.00646) %>% 
  select(X)

write.table(H.punctata_HYPU_WGD_paralogs, quote = F, file = "Hypolepis_punctata_HYPU_WGD_paralogs.tsv")

#### Lepisorus albertii LEAL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.albertii_LEAL <- read.delim("ks_distributions/Lepisorus_albertii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.albertii_LEAL_filt <- L.albertii_LEAL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.albertii_LEAL_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.albertii_LEAL_filt, mapping = aes(x=Ks, ..scaled..*1300), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lepisorus albertii (LEAL)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.albertii_LEAL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.albertii_LEAL_filt.0 <- L.albertii_LEAL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.albertii_LEAL_normalmixEM <- normalmixEM(L.albertii_LEAL_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(L.albertii_LEAL_normalmixEM )

#summary of normalmixEM object:
#        comp 1   comp 2 
#lambda 0.184250 0.815750
#mu     0.251473 2.280188
#sigma  0.201125 0.956878
#loglik at estimate:  -9102.174 

# Generate list of paralogs 

L.albertii_LEAL_WGD_paralogs <-L.albertii_LEAL %>% 
  filter(Ks > 1.32331) %>% 
  filter(Ks < 3.237066) %>% 
  select(X)

write.table(L.albertii_LEAL_WGD_paralogs, quote = F, file = "Lepisorus_albertii_LEAL_WGD_paralogs.tsv")

#### Lepisorus astrolepis LEAS ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.astrolepis_LEAS <- read.delim("ks_distributions/Lepisorus_astrolepis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.astrolepis_LEAS_filt <- L.astrolepis_LEAS %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.astrolepis_LEAS_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.astrolepis_LEAS_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lepisorus astrolepis (LEAS)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.astrolepis_LEAS.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.astrolepis_LEAS_filt.0 <- L.astrolepis_LEAS_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.astrolepis_LEAS_normalmixEM <- normalmixEM(L.astrolepis_LEAS_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak, comp 1 is overfit? 
summary(L.astrolepis_LEAS_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.218879 0.781121
#mu     0.265185 2.262612
#sigma  0.211045 0.952315
#loglik at estimate:  -6968.369

# Generate list of paralogs 

L.astrolepis_LEAS_WGD_paralogs <-L.astrolepis_LEAS %>% 
  filter(Ks > 1.310297) %>% 
  filter(Ks < 3.214927) %>% 
  select(X)

write.table(L.astrolepis_LEAS_WGD_paralogs, quote = F, file = "Lepisorus_astrolepis_LEAS_WGD_paralogs.tsv")

#### Leptochilus cantoniensis LECA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.cantoniensis_LECA <- read.delim("ks_distributions/Leptochilus_cantoniensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.cantoniensis_LECA_filt <- L.cantoniensis_LECA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.cantoniensis_LECA_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.cantoniensis_LECA_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Leptochilus cantoniensis (LECA)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.cantoniensis_LECA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.cantoniensis_LECA_filt.0 <- L.cantoniensis_LECA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.cantoniensis_LECA_normalmixEM <- normalmixEM(L.cantoniensis_LECA_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(L.cantoniensis_LECA_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.888907 0.111093
#mu     2.216157 0.190601
#sigma  0.958760 0.148483
#loglik at estimate:  -9476.81
# Generate list of paralogs 

L.cantoniensis_LECA_WGD_paralogs <-L.cantoniensis_LECA %>% 
  filter(Ks > 1.257397) %>% 
  filter(Ks < 3.174917) %>% 
  select(X)

write.table(L.cantoniensis_LECA_WGD_paralogs, quote = F, file = "Leptochilus_cantoniensis_LECA_WGD_paralogs.tsv")

#### Leptochilus ellipticus LEEL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.ellipticus_LEEL <- read.delim("ks_distributions/Leptochilus_ellipticus.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.ellipticus_LEEL_filt <- L.ellipticus_LEEL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.ellipticus_LEEL_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.ellipticus_LEEL_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Leptochilus ellipticus (LEEL)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.ellipticus_LEEL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.ellipticus_LEEL_filt.0 <- L.ellipticus_LEEL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.ellipticus_LEEL_normalmixEM <- normalmixEM(L.ellipticus_LEEL_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD peak
summary(L.ellipticus_LEEL_normalmixEM )

#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.153628 0.846372
#mu     0.300744 2.286595
#sigma  0.212805 0.921936
#loglik at estimate:  -6698.288

# Generate list of paralogs 

L.ellipticus_LEEL_WGD_paralogs <-L.ellipticus_LEEL %>% 
  filter(Ks > 1.364569) %>% 
  filter(Ks < 3.208531) %>% 
  select(X)

write.table(L.ellipticus_LEEL_WGD_paralogs, quote = F, file = "Leptochilus_ellipticus_LEEL_WGD_paralogs.tsv")

#### Leucostegia immersa LEIM ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.immersa_LEIM <- read.delim("ks_distributions/Leucostegia_immersa_LEIM.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.immersa_LEIM_filt <- L.immersa_LEIM %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.immersa_LEIM_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.immersa_LEIM_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Leucostegia immersa (LEIM)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.immersa_LEIM.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.immersa_LEIM_filt.0 <- L.immersa_LEIM_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.immersa_LEIM_normalmixEM <- normalmixEM(L.immersa_LEIM_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(L.immersa_LEIM_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.103024 0.896976
#mu     0.158555 2.156675
#sigma  0.125867 0.969222
#loglik at estimate:  -6657.821 

# Generate list of paralogs 

L.immersa_LEIM_WGD_paralogs <-L.immersa_LEIM %>% 
  filter(Ks > 1.187453) %>% 
  filter(Ks < 3.125897) %>% 
  select(X)

write.table(L.immersa_LEIM_WGD_paralogs, quote = F, file = "Leucostegia_immersa_LEIM_WGD_paralogs.tsv")

#### Leucostegia immersa WGTU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.immersa_WGTU <- read.delim("ks_distributions/Leucostegia_immersa_WGTU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.immersa_WGTU_filt <- L.immersa_WGTU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.immersa_WGTU_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.immersa_WGTU_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Leucostegia immersa (WGTU)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.immersa_WGTU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.immersa_WGTU_filt.0 <- L.immersa_WGTU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.immersa_WGTU_normalmixEM <- normalmixEM(L.immersa_WGTU_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(L.immersa_WGTU_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.112810 0.887190
#mu     0.124191 2.155607
#sigma  0.109277 0.985063
#loglik at estimate:  -8589.646

# Generate list of paralogs 

L.immersa_WGTU_WGD_paralogs <-L.immersa_WGTU %>% 
  filter(Ks > 1.170544) %>% 
  filter(Ks < 3.14067) %>% 
  select(X)

write.table(L.immersa_WGTU_WGD_paralogs, quote = F, file = "Leucostegia_immersa_WGTU_WGD_paralogs.tsv")

#### Lindsaea heterophylla LIHE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.heterophylla_LIHE <- read.delim("ks_distributions/Lindsaea_heterophylla.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.heterophylla_LIHE_filt <- L.heterophylla_LIHE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.heterophylla_LIHE_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.heterophylla_LIHE_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lindsaea heterophylla (LIHE)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.heterophylla_LIHE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.heterophylla_LIHE_filt.0 <- L.heterophylla_LIHE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.heterophylla_LIHE_normalmixEM <- normalmixEM(L.heterophylla_LIHE_filt.0$Ks, k=2)
# Fits to 2 WGD peaks
summary(L.heterophylla_LIHE_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.422747 0.577253
#mu     0.720338 2.556214
#sigma  0.425086 0.790241
#loglik at estimate:  -8493.723

# Generate list of paralogs 

L.heterophylla_LIHE_WGD_paralogs_peak1 <-L.heterophylla_LIHE %>% 
  filter(Ks > 0.295252) %>% 
  filter(Ks < 1.145424) %>% 
  select(X)

write.table(L.heterophylla_LIHE_WGD_paralogs_peak1, quote = F, file = "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak1.tsv")

L.heterophylla_LIHE_WGD_paralogs_peak2 <-L.heterophylla_LIHE %>% 
  filter(Ks > 1.765973) %>% 
  filter(Ks < 3.346455) %>% 
  select(X)

write.table(L.heterophylla_LIHE_WGD_paralogs_peak2, quote = F, file = "Lindsaea_heterophylla_LIHE_WGD_paralogs_peak2.tsv")

#### Lindsaea linearis NOKI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.linearis_NOKI <- read.delim("ks_distributions/Lindsaea_linearis_NOKI.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.linearis_NOKI_filt <- L.linearis_NOKI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.linearis_NOKI_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.linearis_NOKI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lindsaea linearis (NOKI)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.linearis_NOKI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.linearis_NOKI_filt.0 <- L.linearis_NOKI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.linearis_NOKI_normalmixEM <- normalmixEM(L.linearis_NOKI_filt.0$Ks, k=2)
# Fits to 2 WGD peaks- first component is a WGD, not recent dups 
summary(L.linearis_NOKI_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.439908 0.560092
#mu     0.602108 2.425540
#sigma  0.411539 0.846740
#loglik at estimate:  -9497.57

# Generate list of paralogs 

L.linearis_NOKI_WGD_paralogs_peak1 <-L.linearis_NOKI %>% 
  filter(Ks > 0.190569) %>% 
  filter(Ks < 1.013647) %>% 
  select(X)

write.table(L.linearis_NOKI_WGD_paralogs_peak1, quote = F, file = "Lindsaea_linearis_NOKI_WGD_paralogs_peak1.tsv")

L.linearis_NOKI_WGD_paralogs_peak2 <-L.linearis_NOKI %>% 
  filter(Ks > 1.5788) %>% 
  filter(Ks < 3.272228) %>% 
  select(X)

write.table(L.linearis_NOKI_WGD_paralogs_peak2, quote = F, file = "Lindsaea_linearis_NOKI_WGD_paralogs_peak2.tsv")

#### Lindsaea microphylla YIXP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.microphylla_YIXP <- read.delim("ks_distributions/Lindsaea_microphylla_YIXP.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.microphylla_YIXP_filt <- L.microphylla_YIXP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.microphylla_YIXP_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.microphylla_YIXP_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lindsaea microphylla (YIXP)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.microphylla_YIXP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.microphylla_YIXP_filt.0 <- L.microphylla_YIXP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.microphylla_YIXP_normalmixEM <- normalmixEM(L.microphylla_YIXP_filt.0$Ks, k=2)
# Fits to 2 WGD peaks- first component is a WGD, not recent dups 
summary(L.microphylla_YIXP_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.363997 0.636003
#mu     0.700350 2.548843
#sigma  0.464089 0.793773
#loglik at estimate:  -8408.116

# Generate list of paralogs 

L.microphylla_YIXP_WGD_paralogs_peak1 <-L.microphylla_YIXP %>% 
  filter(Ks > 0.236261) %>% 
  filter(Ks < 1.164439) %>% 
  select(X)

write.table(L.microphylla_YIXP_WGD_paralogs_peak1, quote = F, file = "Lindsaea_microphylla_YIXP_WGD_paralogs_peak1.tsv")

L.microphylla_YIXP_WGD_paralogs_peak2 <-L.microphylla_YIXP %>% 
  filter(Ks > 1.75507) %>% 
  filter(Ks < 3.342616) %>% 
  select(X)

write.table(L.microphylla_YIXP_WGD_paralogs_peak2, quote = F, file = "Lindsaea_microphylla_YIXP_WGD_paralogs_peak2.tsv")

#### Lomagramma matthewii LOMA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.matthewii_LOMA <- read.delim("ks_distributions/Lomagramma_matthewii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.matthewii_LOMA_filt <- L.matthewii_LOMA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.matthewii_LOMA_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.matthewii_LOMA_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lomagramma matthewii (LOMA)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.matthewii_LOMA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.matthewii_LOMA_filt.0 <- L.matthewii_LOMA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.matthewii_LOMA_normalmixEM <- normalmixEM(L.matthewii_LOMA_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(L.matthewii_LOMA_normalmixEM )
#summary of normalmixEM object:
#          comp 1   comp 2
#lambda 0.116770 0.883230
#mu     0.147294 2.218476
#sigma  0.119906 0.960661
#loglik at estimate:  -6853.146

# Generate list of paralogs 

L.matthewii_LOMA_WGD_paralogs <-L.matthewii_LOMA %>% 
  filter(Ks > 1.257815) %>% 
  filter(Ks < 3.179137) %>% 
  select(X)

write.table(L.matthewii_LOMA_WGD_paralogs, quote = F, file = "Lomagramma_matthewii_LOMA_WGD_paralogs.tsv")

#### Lomagramma sumatrana LOSU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.sumatrana_LOSU <- read.delim("ks_distributions/Lomagramma_sumatrana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.sumatrana_LOSU_filt <- L.sumatrana_LOSU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.sumatrana_LOSU_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.sumatrana_LOSU_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lomagramma sumatrana (LOSU)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.sumatrana_LOSU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.sumatrana_LOSU_filt.0 <- L.sumatrana_LOSU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.sumatrana_LOSU_normalmixEM <- normalmixEM(L.sumatrana_LOSU_filt.0$Ks, k=2)
# Fits to 1 WGD peaks
summary(L.sumatrana_LOSU_normalmixEM )
#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.115592 0.884408
#mu     0.177274 2.204308
#sigma  0.154183 0.964503
#loglik at estimate:  -6803.563 

# Generate list of paralogs 

L.sumatrana_LOSU_WGD_paralogs <-L.sumatrana_LOSU %>% 
  filter(Ks > 1.239805) %>% 
  filter(Ks < 3.168811) %>% 
  select(X)

write.table(L.sumatrana_LOSU_WGD_paralogs, quote = F, file = "Lomagramma_sumatrana_LOSU_WGD_paralogs.tsv")

#### Lomariopsis boninensis LOBO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.boninensis_LOBO <- read.delim("ks_distributions/Lomariopsis_boninensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.boninensis_LOBO_filt <- L.boninensis_LOBO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.boninensis_LOBO_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.boninensis_LOBO_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lomariopsis boninensis (LOBO)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.boninensis_LOBO.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD
 
#### Lomariopsis spectabilis LOSP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.spectabilis_LOSP <- read.delim("ks_distributions/Lomariopsis_spectabilis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.spectabilis_LOSP_filt <- L.spectabilis_LOSP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.spectabilis_LOSP_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.spectabilis_LOSP_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lomariopsis spectabilis (LOSP)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.spectabilis_LOSP.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Lonchitis hirsuta VVRN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.hirsuta_VVRN <- read.delim("ks_distributions/Lonchitis_hirsuta_VVRN.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.hirsuta_VVRN_filt <- L.hirsuta_VVRN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.hirsuta_VVRN_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.hirsuta_VVRN_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lonchitis hirsuta (VVRN)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.hirsuta_VVRN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.hirsuta_VVRN_filt.0 <- L.hirsuta_VVRN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.hirsuta_VVRN_normalmixEM <- normalmixEM(L.hirsuta_VVRN_filt.0$Ks, k=3)
# Fits to 2 WGD peaks- ignore last component 
summary(L.hirsuta_VVRN_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.312070 0.568234 0.119696
#mu     0.353514 1.920793 3.478851
#sigma  0.239761 0.782836 0.308298
#loglik at estimate:  -7950.963

# Generate list of paralogs 

L.hirsuta_VVRN_WGD_paralogs_peak1 <-L.hirsuta_VVRN %>% 
  filter(Ks > 0.113753) %>% 
  filter(Ks < 0.593275) %>% 
  select(X)

write.table(L.hirsuta_VVRN_WGD_paralogs_peak1, quote = F, file = "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak1.tsv")

L.hirsuta_VVRN_WGD_paralogs_peak2 <-L.hirsuta_VVRN %>% 
  filter(Ks > 1.137957) %>% 
  filter(Ks < 2.703629) %>% 
  select(X)

write.table(L.hirsuta_VVRN_WGD_paralogs_peak2, quote = F, file = "Lonchitis_hirsuta_VVRN_WGD_paralogs_peak2.tsv")

#### Loxogramma biformis LOBI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.biformis_LOBI <- read.delim("ks_distributions/Loxogramme_biformis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.biformis_LOBI_filt <- L.biformis_LOBI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.biformis_LOBI_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.biformis_LOBI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Loxogramme biformis (LOBI)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.biformis_LOBI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.biformis_LOBI_filt.0 <- L.biformis_LOBI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.biformis_LOBI_normalmixEM <- normalmixEM(L.biformis_LOBI_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(L.biformis_LOBI_normalmixEM ) 
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.143612 0.856388
#mu     0.254709 2.340844
#sigma  0.201976 0.918548
#loglik at estimate:  -6683.08 

# Generate list of paralogs 

L.biformis_LOBI_WGD_paralogs <- L.biformis_LOBI %>% 
  filter(Ks > 1.422296) %>% 
  filter(Ks < 3.259392) %>% 
  select(X)

write.table(L.biformis_LOBI_WGD_paralogs, quote = F, file = "Loxogramme_biformis_LOBI_WGD_paralogs.tsv")

#### Loxogramme chinensis LOCH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.chinensis_LOCH <- read.delim("ks_distributions/Loxogramme_chinensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.chinensis_LOCH_filt <- L.chinensis_LOCH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.chinensis_LOCH_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.chinensis_LOCH_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Loxogramme chinensis (LOCH)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.chinensis_LOCH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.chinensis_LOCH_filt.0 <- L.chinensis_LOCH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.chinensis_LOCH_normalmixEM <- normalmixEM(L.chinensis_LOCH_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(L.chinensis_LOCH_normalmixEM ) 
#summary of normalmixEM object:
#  comp 1  comp 2
#lambda 0.182530 0.81747
#mu     0.131606 2.12463
#sigma  0.105177 1.01878
#loglik at estimate:  -6519.629

# Generate list of paralogs 

L.chinensis_LOCH_WGD_paralogs <- L.chinensis_LOCH %>% 
  filter(Ks > 1.10585) %>% 
  filter(Ks < 3.14341) %>% 
  select(X)

write.table(L.chinensis_LOCH_WGD_paralogs, quote = F, file = "Loxogramme_chinensis_LOCH_WGD_paralogs.tsv")

#### Lygodium flexuosum LYFL####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.flexuosum_LYFL <- read.delim("ks_distributions/Lygodium_flexuosum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.flexuosum_LYFL_filt <- L.flexuosum_LYFL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.flexuosum_LYFL_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.flexuosum_LYFL_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lygodium flexuosum (LYFL)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.flexuosum_LYFL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.flexuosum_LYFL_filt.0 <- L.flexuosum_LYFL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.flexuosum_LYFL_normalmixEM <- normalmixEM(L.flexuosum_LYFL_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(L.flexuosum_LYFL_normalmixEM ) 
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.182745 0.817255
#mu     0.129411 2.140283
#sigma  0.120438 1.002557
#loglik at estimate:  -6319.34 

# Generate list of paralogs 

L.flexuosum_LYFL_WGD_paralogs <- L.flexuosum_LYFL %>% 
  filter(Ks > 1.137726) %>% 
  filter(Ks < 3.14284) %>% 
  select(X)

write.table(L.flexuosum_LYFL_WGD_paralogs, quote = F, file = "Lygodium_flexuosum_LYFL_WGD_paralogs.tsv")

#### Lygodium japonicum LYJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.japonicum_LYJA <- read.delim("ks_distributions/Lygodium_japonicum_LYJA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.japonicum_LYJA_filt <- L.japonicum_LYJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.japonicum_LYJA_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.japonicum_LYJA_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lygodium japonicum (LYJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.japonicum_LYJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.japonicum_LYJA_filt.0 <- L.japonicum_LYJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.japonicum_LYJA_normalmixEM <- normalmixEM(L.japonicum_LYJA_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(L.japonicum_LYJA_normalmixEM ) 
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0966811 0.903319
#mu     0.1272349 2.133877
#sigma  0.0990098 0.984444
#loglik at estimate:  -7570.03 

# Generate list of paralogs 

L.japonicum_LYJA_WGD_paralogs <- L.japonicum_LYJA %>% 
  filter(Ks > 1.149433) %>% 
  filter(Ks < 3.118321) %>% 
  select(X)

write.table(L.japonicum_LYJA_WGD_paralogs, quote = F, file = "Lygodium_japonicum_LYJA_WGD_paralogs.tsv")

#### Lygodium japonicum PBUU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

L.japonicum_PBUU <- read.delim("ks_distributions/Lygodium_japonicum_PBUU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

L.japonicum_PBUU_filt <- L.japonicum_PBUU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(L.japonicum_PBUU_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(L.japonicum_PBUU_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Lygodium japoncium (PBUU)") + theme(plot.title = element_text(face = "italic"))

ggsave("L.japonicum_PBUU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

L.japonicum_PBUU_filt.0 <- L.japonicum_PBUU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

L.japonicum_PBUU_normalmixEM <- normalmixEM(L.japonicum_PBUU_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(L.japonicum_PBUU_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.2055860 0.794414
#mu     0.0933235 2.190893
#sigma  0.0915093 0.994686
#loglik at estimate:  -4125.392 

# Generate list of paralogs 

L.japonicum_PBUU_WGD_paralogs <- L.japonicum_PBUU %>% 
  filter(Ks > 1.196207) %>% 
  filter(Ks < 3.185579) %>% 
  select(X)

write.table(L.japonicum_PBUU_WGD_paralogs, quote = F, file = "Lygodium_japonicum_PBUU_WGD_paralogs.tsv")

#### Marsilea quadrifolia MAQF ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.quadrifolia_MAQF <- read.delim("ks_distributions/Marsilea_quadrifolia_MAQF.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.quadrifolia_MAQF_filt <- M.quadrifolia_MAQF %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.quadrifolia_MAQF_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.quadrifolia_MAQF_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Marsilea quadrifolia (MAQF)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.quadrifolia_MAQF.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.quadrifolia_MAQF_filt.0 <- M.quadrifolia_MAQF_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.quadrifolia_MAQF_normalmixEM <- normalmixEM(M.quadrifolia_MAQF_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(M.quadrifolia_MAQF_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1756854 0.824315
#mu     0.0834925 2.131091
#sigma  0.0760720 1.051117
#loglik at estimate:  -5744.668 

# Generate list of paralogs 

M.quadrifolia_MAQF_WGD_paralogs <- M.quadrifolia_MAQF %>% 
  filter(Ks > 1.079974) %>% 
  filter(Ks < 3.182208) %>% 
  select(X)

write.table(M.quadrifolia_MAQF_WGD_paralogs, quote = F, file = "Marsilea_quadrifolia_MAQF_WGD_paralogs.tsv")

#### Marsilea quadrifolia MAQU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.quadrifolia_MAQU <- read.delim("ks_distributions/Marsilea_quadrifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.quadrifolia_MAQU_filt <- M.quadrifolia_MAQU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.quadrifolia_MAQU_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.quadrifolia_MAQU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Marsilea quadrifolia (MAQU)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.quadrifolia_MAQU.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Matteuccia struthiopteris MAST ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.struthiopteris_MAST <- read.delim("ks_distributions/Matteuccia_struthiopteris.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.struthiopteris_MAST_filt <- M.struthiopteris_MAST %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.struthiopteris_MAST_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.struthiopteris_MAST_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Matteuccia struthiopteris (MAST)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.struthiopteris_MAST.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.struthiopteris_MAST_filt.0 <- M.struthiopteris_MAST_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.struthiopteris_MAST_normalmixEM <- normalmixEM(M.struthiopteris_MAST_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(M.struthiopteris_MAST_normalmixEM )
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.118131 0.881869
#mu     0.228399 2.118557
#sigma  0.154281 0.952301
#loglik at estimate:  -8332.283 

# Generate list of paralogs 

M.struthiopteris_MAST_WGD_paralogs <- M.struthiopteris_MAST %>% 
  filter(Ks > 1.166256) %>% 
  filter(Ks < 3.070858) %>% 
  select(X)

write.table(M.struthiopteris_MAST_WGD_paralogs, quote = F, file = "Matteuccia_struthiopteris_MAST_WGD_paralogs.tsv")

#### Microlepia hookeriana MIHO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.hookeriana_MIHO <- read.delim("ks_distributions/Microlepia_hookeriana.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.hookeriana_MIHO_filt <- M.hookeriana_MIHO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.hookeriana_MIHO_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.hookeriana_MIHO_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Microlepia hookeriana (MIHO)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.hookeriana_MIHO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

# No evidence of WGD 

#### Microlepia marginata MIMA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.marginata_MIMA <- read.delim("ks_distributions/Microlepia_marginata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.marginata_MIMA_filt <- M.marginata_MIMA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.marginata_MIMA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.marginata_MIMA_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Microlepia marginata (MIMA)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.marginata_MIMA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.marginata_MIMA_filt.0 <- M.marginata_MIMA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.marginata_MIMA_normalmixEM <- normalmixEM(M.marginata_MIMA_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(M.marginata_MIMA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.822370 0.177630
#mu     1.518186 3.343603
#sigma  0.913802 0.403995
#loglik at estimate:  -10603.2  

# Generate list of paralogs 

M.marginata_MIMA_WGD_paralogs <- M.marginata_MIMA %>% 
  filter(Ks > 0.604384) %>% 
  filter(Ks < 2.431988) %>% 
  select(X)

write.table(M.marginata_MIMA_WGD_paralogs, quote = F, file = "Microlepia_marginata_MIMA_WGD_paralogs.tsv")

#### Microlepia platyphylla MIPL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.platyphylla_MIPL <- read.delim("ks_distributions/Microlepia_platyphylla.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.platyphylla_MIPL_filt <- M.platyphylla_MIPL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.platyphylla_MIPL_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.platyphylla_MIPL_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Microlepia platyphylla (MIPL)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.platyphylla_MIPL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.platyphylla_MIPL_filt.0 <- M.platyphylla_MIPL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.platyphylla_MIPL_normalmixEM <- normalmixEM(M.platyphylla_MIPL_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(M.platyphylla_MIPL_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.113913 0.886087
#mu     0.181753 2.192921
#sigma  0.147443 0.952081
#loglik at estimate:  -8470.022 

# Generate list of paralogs 

M.platyphylla_MIPL_WGD_paralogs <- M.platyphylla_MIPL %>% 
  filter(Ks > 1.24084) %>% 
  filter(Ks < 3.145002) %>% 
  select(X)

write.table(M.platyphylla_MIPL_WGD_paralogs, quote = F, file = "Microlepia_platyphylla_MIPL_WGD_paralogs.tsv")

#### Microlepia speluncae MISP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.speluncae_MISP<- read.delim("ks_distributions/Microlepia_speluncae.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.speluncae_MISP_filt <- M.speluncae_MISP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.speluncae_MISP_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.speluncae_MISP_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Microlepia speluncae (MISP)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.speluncae_MISP.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Microsorum scolopendria MISC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.scolopendria_MISC <- read.delim("ks_distributions/Microsorum_scolopendria.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.scolopendria_MISC_filt <- M.scolopendria_MISC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.scolopendria_MISC_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.scolopendria_MISC_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Microsorum scolopendria (MISC)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.scolopendria_MISC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.scolopendria_MISC_filt.0 <- M.scolopendria_MISC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.scolopendria_MISC_normalmixEM <- normalmixEM(M.scolopendria_MISC_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(M.scolopendria_MISC_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1419979 0.858002
#mu     0.1052764 1.989689
#sigma  0.0857082 1.034198
#loglik at estimate:  -9880.568

# Generate list of paralogs 

M.scolopendria_MISC_WGD_paralogs <- M.scolopendria_MISC %>% 
  filter(Ks > 0.955491) %>% 
  filter(Ks < 3.023887) %>% 
  select(X)

write.table(M.scolopendria_MISC_WGD_paralogs, quote = F, file = "Microsorum_scolopendria_MISC_WGD_paralogs.tsv")

#### Monachosorum flagellare MOFL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.flagellare_MOFL <- read.delim("ks_distributions/Monachosorum_flagellare.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.flagellare_MOFL_filt <- M.flagellare_MOFL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.flagellare_MOFL_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.flagellare_MOFL_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Monachosorum flagellare (MOFL)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.flagellare_MOFL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.flagellare_MOFL_filt.0 <- M.flagellare_MOFL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.flagellare_MOFL_normalmixEM <- normalmixEM(M.flagellare_MOFL_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(M.flagellare_MOFL_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.148809 0.851191
#mu     0.108367 2.016109
#sigma  0.093410 0.978949
#loglik at estimate:  -5818.22

# Generate list of paralogs 

M.flagellare_MOFL_WGD_paralogs <- M.flagellare_MOFL %>% 
  filter(Ks > 1.03716) %>% 
  filter(Ks < 2.995058) %>% 
  select(X)

write.table(M.flagellare_MOFL_WGD_paralogs, quote = F, file = "Monachosorum_flagellare_MOFL_WGD_paralogs.tsv")

#### Monachosorum henryi MOHE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.henryi_MOHE <- read.delim("ks_distributions/Monachosorum_henryi.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.henryi_MOHE_filt <- M.henryi_MOHE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.henryi_MOHE_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.henryi_MOHE_filt, mapping = aes(x=Ks, ..scaled..*1300), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Monachosorum henryi (MOHE)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.henryi_MOHE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.henryi_MOHE_filt.0 <- M.henryi_MOHE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.henryi_MOHE_normalmixEM <- normalmixEM(M.henryi_MOHE_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(M.henryi_MOHE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1651581 0.834842
#mu     0.1035154 2.049784
#sigma  0.0843835 1.018901
#loglik at estimate:  -8215.263

# Generate list of paralogs 

M.henryi_MOHE_WGD_paralogs <- M.henryi_MOHE %>% 
  filter(Ks > 1.030883) %>% 
  filter(Ks < 3.068685) %>% 
  select(X)

write.table(M.henryi_MOHE_WGD_paralogs, quote = F, file = "Monachosorum_henryi_MOHE_WGD_paralogs.tsv")

#### Monachosorum maximowiczii MOMA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.maximowiczii_MOMA <- read.delim("ks_distributions/Monachosorum_maximowiczii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.maximowiczii_MOMA_filt <- M.maximowiczii_MOMA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.maximowiczii_MOMA_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.maximowiczii_MOMA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Monachosorum maximowiczii (MOMA)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.maximowiczii_MOMA.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Myriopteris rufa GSXD ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

M.rufa_GSXD <- read.delim("ks_distributions/Myriopteris_rufa_GSXD.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

M.rufa_GSXD_filt <- M.rufa_GSXD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(M.rufa_GSXD_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(M.rufa_GSXD_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Myriopteris rufa (GSXD)") + theme(plot.title = element_text(face = "italic"))

ggsave("M.rufa_GSXD.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

M.rufa_GSXD_filt.0 <- M.rufa_GSXD_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

M.rufa_GSXD_normalmixEM <- normalmixEM(M.rufa_GSXD_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(M.rufa_GSXD_normalmixEM )
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1352497 0.864750
#mu     0.1241726 2.216651
#sigma  0.0948084 0.988503
#loglik at estimate:  -6111.004 

# Generate list of paralogs 

M.rufa_GSXD_WGD_paralogs <- M.rufa_GSXD %>% 
  filter(Ks > 1.228148) %>% 
  filter(Ks < 3.205154) %>% 
  select(X)

write.table(M.rufa_GSXD_WGD_paralogs, quote = F, file = "Myriopteris_rufa_GSXD_WGD_paralogs.tsv")

#### Nephrolepis biserrata NEBI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

N.biserrata_NEBI <- read.delim("ks_distributions/Nephrolepis_biserrata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

N.biserrata_NEBI_filt <- N.biserrata_NEBI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(N.biserrata_NEBI_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(N.biserrata_NEBI_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Nephrolepis biserrata (NEBI)") + theme(plot.title = element_text(face = "italic"))

ggsave("N.biserrata_NEBI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

N.biserrata_NEBI_filt.0 <- N.biserrata_NEBI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

N.biserrata_NEBI_normalmixEM <- normalmixEM(N.biserrata_NEBI_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(N.biserrata_NEBI_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0858747 0.914125
#mu     0.0961484 2.222618
#sigma  0.0916671 0.962967
#loglik at estimate:  -6616.198 

# Generate list of paralogs 

N.biserrata_NEBI_WGD_paralogs <- N.biserrata_NEBI %>% 
  filter(Ks > 1.259651) %>% 
  filter(Ks < 3.185585) %>% 
  select(X)

write.table(N.biserrata_NEBI_WGD_paralogs, quote = F, file = "Nephrolepis_biserrata_NEBI_WGD_paralogs.tsv")

#### Nephrolepis cordifolia NECO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

N.cordifolia_NECO <- read.delim("ks_distributions/Nephrolepis_cordifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

N.cordifolia_NECO_filt <- N.cordifolia_NECO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(N.cordifolia_NECO_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(N.cordifolia_NECO_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Nephrolepis cordifolia (NECO)") + theme(plot.title = element_text(face = "italic"))

ggsave("N.cordifolia_NECO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

N.cordifolia_NECO_filt.0 <- N.cordifolia_NECO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

N.cordifolia_NECO_normalmixEM <- normalmixEM(N.cordifolia_NECO_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(N.cordifolia_NECO_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.899357 0.100643
#mu     2.234918 0.200559
#sigma  0.944579 0.146843
#loglik at estimate:  -6635.071 

# Generate list of paralogs 

N.cordifolia_NECO_WGD_paralogs <- N.cordifolia_NECO %>% 
  filter(Ks > 1.290339) %>% 
  filter(Ks < 3.179497) %>% 
  select(X)

write.table(N.cordifolia_NECO_WGD_paralogs, quote = F, file = "Nephrolepis_cordifolia_NECO_WGD_paralogs.tsv")

#### Notholaena montieliae YCKE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

N.montieliae_YCKE <- read.delim("ks_distributions/Notholaena_montieliae_YCKE.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

N.montieliae_YCKE_filt <- N.montieliae_YCKE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(N.montieliae_YCKE_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(N.montieliae_YCKE_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Notholaena montieliae (YCKE)") + theme(plot.title = element_text(face = "italic"))

ggsave("N.montieliae_YCKE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

N.montieliae_YCKE_filt.0 <- N.montieliae_YCKE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

N.montieliae_YCKE_normalmixEM <- normalmixEM(N.montieliae_YCKE_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(N.montieliae_YCKE_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.123626 0.876374
#mu     0.123343 2.164982
#sigma  0.106375 0.999527
#loglik at estimate:  -5327.107

# Generate list of paralogs 

N.montieliae_YCKE_WGD_paralogs <- N.montieliae_YCKE %>% 
  filter(Ks > 1.165455) %>% 
  filter(Ks < 3.164509) %>% 
  select(X)

write.table(N.montieliae_YCKE_WGD_paralogs, quote = F, file = "Notholaena_montieliae_YCKE_WGD_paralogs.tsv")

#### Odontosoria chinensis ODCH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.chinensis_ODCH <- read.delim("ks_distributions/Odontosira_chinensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.chinensis_ODCH_filt <- O.chinensis_ODCH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.chinensis_ODCH_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.chinensis_ODCH_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Odontosoria chinensis (ODCH)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.chinensis_ODCH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.chinensis_ODCH_filt.0 <- O.chinensis_ODCH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.chinensis_ODCH_normalmixEM <- normalmixEM(O.chinensis_ODCH_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(O.chinensis_ODCH_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.408014 0.591986
#mu     0.603334 2.321527
#sigma  0.389831 0.874429
#loglik at estimate:  -11107.16 

# Generate list of paralogs 

O.chinensis_ODCH_WGD_paralogs <- O.chinensis_ODCH %>% 
  filter(Ks > 0.213503) %>% 
  filter(Ks < 0.993165) %>% 
  select(X)

write.table(O.chinensis_ODCH_WGD_paralogs, quote = F, file = "Odontosira_chinensis_ODCH_WGD_paralogs.tsv")

#### Oleandra hainanensis OLHA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.hainanensis_OLHA <- read.delim("ks_distributions/Oleandra_hainanensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.hainanensis_OLHA_filt <- O.hainanensis_OLHA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.hainanensis_OLHA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.hainanensis_OLHA_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Oleandra hainanensis (OLHA)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.hainanensis_OLHA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.hainanensis_OLHA_filt.0 <- O.hainanensis_OLHA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.hainanensis_OLHA_normalmixEM <- normalmixEM(O.hainanensis_OLHA_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(O.hainanensis_OLHA_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.831541 0.168459
#mu     1.760283 3.444440
#sigma  0.946831 0.343009
#loglik at estimate:  -7256.099 

# Generate list of paralogs 

O.hainanensis_OLHA_WGD_paralogs <- O.hainanensis_OLHA %>% 
  filter(Ks > 0.813452) %>% 
  filter(Ks < 2.707114) %>% 
  select(X)

write.table(O.hainanensis_OLHA_WGD_paralogs, quote = F, file = "Oleandra_hainanensis_OLHA_WGD_paralogs.tsv")

#### Oleandra musifolia OLMU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.musifolia_OLMU <- read.delim("ks_distributions/Oleandra_musifolia.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.musifolia_OLMU_filt <- O.musifolia_OLMU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.musifolia_OLMU_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.musifolia_OLMU_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Oleandra musifolia (OLMU)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.musifolia_OLMU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.musifolia_OLMU_filt.0 <- O.musifolia_OLMU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.musifolia_OLMU_normalmixEM <- normalmixEM(O.musifolia_OLMU_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(O.musifolia_OLMU_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.141500 0.858500
#mu     0.229993 2.196089
#sigma  0.151825 0.940295
#loglik at estimate:  -8128.545 

# Generate list of paralogs 

O.musifolia_OLMU_WGD_paralogs <- O.musifolia_OLMU %>% 
  filter(Ks > 1.255794) %>% 
  filter(Ks < 3.136384) %>% 
  select(X)

write.table(O.musifolia_OLMU_WGD_paralogs, quote = F, file = "Oleandra_musifolia_OLMU_WGD_paralogs.tsv")

#### Oleandra sp. OLSP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.sp_OLSP <- read.delim("ks_distributions/Oleandra_sp._XQ-2018.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.sp_OLSP_filt <- O.sp_OLSP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.sp_OLSP_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.sp_OLSP_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Oleandra sp. (OLSP)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.sp_OLSP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.sp_OLSP_filt.0 <- O.sp_OLSP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.sp_OLSP_normalmixEM <- normalmixEM(O.sp_OLSP_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(O.sp_OLSP_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.119406 0.880594
#mu     0.181434 2.184436
#sigma  0.139260 0.953678
#loglik at estimate:  -8067.961

# Generate list of paralogs 

O.sp_OLSP_WGD_paralogs <- O.sp_OLSP %>% 
  filter(Ks > 1.230758) %>% 
  filter(Ks < 3.138114) %>% 
  select(X)

write.table(O.sp_OLSP_WGD_paralogs, quote = F, file = "Oleandra_sp._XQ-2018_OLSP_WGD_paralogs.tsv")

#### Onoclea sensibilis HTFH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.sensibilis_HTFH <- read.delim("ks_distributions/Onoclea_sensibilis_HTFH.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.sensibilis_HTFH_filt <- O.sensibilis_HTFH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.sensibilis_HTFH_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.sensibilis_HTFH_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Onoclea sensibilis (HTFH)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.sensibilis_HTFH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.sensibilis_HTFH_filt.0 <- O.sensibilis_HTFH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.sensibilis_HTFH_normalmixEM <- normalmixEM(O.sensibilis_HTFH_filt.0$Ks, k=3)
# Fits to 1 WGD peak, ignore last component 
summary(O.sensibilis_HTFH_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.130048 0.728232 0.141719
#mu     0.174408 1.855737 3.492725
#sigma  0.127979 0.796221 0.302341
#loglik at estimate:  -3888.947 

# Generate list of paralogs 

O.sensibilis_HTFH_WGD_paralogs <- O.sensibilis_HTFH  %>% 
  filter(Ks > 1.120843) %>% 
  filter(Ks < 3.110537) %>% 
  select(X)

write.table(O.sensibilis_HTFH_WGD_paralogs, quote = F, file = "Onoclea_sensibilis_HTFH_WGD_paralogs.tsv")

#### Onoclea sensibilis ONSE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.sensibilis_ONSE <- read.delim("ks_distributions/Onoclea_sensibilis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.sensibilis_ONSE_filt <- O.sensibilis_ONSE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.sensibilis_ONSE_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.sensibilis_ONSE_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Onoclea sensibilis (ONSE)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.sensibilis_ONSE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.sensibilis_ONSE_filt.0 <- O.sensibilis_ONSE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.sensibilis_ONSE_normalmixEM <- normalmixEM(O.sensibilis_ONSE_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 is overfit? 
summary(O.sensibilis_ONSE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.112461 0.887539
#mu     0.171244 2.115690
#sigma  0.125470 0.994847
#loglik at estimate:  -9299.885 

# Generate list of paralogs 

O.sensibilis_ONSE_WGD_paralogs <- O.sensibilis_ONSE  %>% 
  filter(Ks > 1.059516) %>% 
  filter(Ks < 2.651958) %>% 
  select(X)

write.table(O.sensibilis_ONSE_WGD_paralogs, quote = F, file = "Onoclea_sensibilis_ONSE_WGD_paralogs.tsv")

#### Onychium japoncium ONJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.japonicum_ONJA <- read.delim("ks_distributions/Onychium_japonicum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.japonicum_ONJA_filt <- O.japonicum_ONJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.japonicum_ONJA_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.japonicum_ONJA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Onychium japonicum (ONJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.japonicum_ONJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.japonicum_ONJA_filt.0 <- O.japonicum_ONJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.japonicum_ONJA_normalmixEM <- normalmixEM(O.japonicum_ONJA_filt.0$Ks, k=2)
# Fits to 1 WGD peak, 
summary(O.japonicum_ONJA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1179192 0.882081
#mu     0.1184120 2.205699
#sigma  0.0991368 0.990761
#loglik at estimate:  -7330.41  

# Generate list of paralogs 

O.japonicum_ONJA_WGD_paralogs <- O.japonicum_ONJA  %>% 
  filter(Ks > 1.214938) %>% 
  filter(Ks < 3.19646) %>% 
  select(X)

write.table(O.japonicum_ONJA_WGD_paralogs, quote = F, file = "Onychium_japonicum_ONJA_WGD_paralogs.tsv")

#### Ophioderma pendula OPPE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.pendula_OPPE <- read.delim("ks_distributions/Ophioderma_pendula.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.pendula_OPPE_filt <- O.pendula_OPPE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.pendula_OPPE_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.pendula_OPPE_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ophioderma pendula (OPPE)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.pendula_OPPE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.pendula_OPPE_filt.0 <- O.pendula_OPPE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.pendula_OPPE_normalmixEM <- normalmixEM(O.pendula_OPPE_filt.0$Ks, k=2)
# Fits to 1 WGD peak,comp 2 is overfit? 
summary(O.pendula_OPPE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.339441 0.660559
#mu     0.387023 2.046387
#sigma  0.222605 0.964960
#loglik at estimate:  -5225.368 

# Generate list of paralogs 

O.pendula_OPPE_WGD_paralogs <- O.pendula_OPPE  %>% 
  filter(Ks > 0.164418) %>% 
  filter(Ks < 0.609628) %>% 
  select(X)

write.table(O.pendula_OPPE_WGD_paralogs, quote = F, file = "Ophioderma_pendula_OPPE_WGD_paralogs.tsv")

#### Ophioglossum thermale OPTH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.thermale_OPTH <- read.delim("ks_distributions/Ophioglossum_thermale.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.thermale_OPTH_filt <- O.thermale_OPTH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.thermale_OPTH_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.thermale_OPTH_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ophioglossum thermale (OPTH)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.thermale_OPTH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.thermale_OPTH_filt.0 <- O.thermale_OPTH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.thermale_OPTH_normalmixEM <- normalmixEM(O.thermale_OPTH_filt.0$Ks, k=3)
# Fits to 1 WGD peak,comp 3 is overfit? 
summary(O.thermale_OPTH_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.269984 0.335249 0.394767
#mu     0.114780 0.894657 2.547129
#sigma  0.085700 0.440976 0.780622
#loglik at estimate:  -6481.53 

# Generate list of paralogs 

O.thermale_OPTH_WGD_paralogs <- O.thermale_OPTH  %>% 
  filter(Ks > 0.453681) %>% 
  filter(Ks < 1.335633) %>% 
  select(X)

write.table(O.thermale_OPTH_WGD_paralogs, quote = F, file = "Ophioglossum_thermale_OPTH_WGD_paralogs.tsv")

#### Ophioglossum vulgatum OPVU ####

O.vulgatum_OPVU <- read.delim("ks_distributions/Ophioglossum_vulgatum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.vulgatum_OPVU_filt <- O.vulgatum_OPVU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.vulgatum_OPVU_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.vulgatum_OPVU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ophioglossum vulgatum (OPVU)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.vulgatum_OPVU.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 
#### Oreogrammitis congener ORCO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.grammitis_ORCO <- read.delim("ks_distributions/Oreogrammitis_congener.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.grammitis_ORCO_filt <- O.grammitis_ORCO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.grammitis_ORCO_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.grammitis_ORCO_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Oreogrammitis congener (ORCO)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.congener_ORCO.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD? 

#### Osmolindsaea ordata OSOD ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.ordata_OSOD <- read.delim("ks_distributions/Osmolindsaea_odorata_OSOD.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.ordata_OSOD_filt <- O.ordata_OSOD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.ordata_OSOD_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.ordata_OSOD_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Osmolindsaea ordata (OSOD)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.ordata_OSOD.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.ordata_OSOD_filt.0 <- O.ordata_OSOD_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.ordata_OSOD_normalmixEM <- normalmixEM(O.ordata_OSOD_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2 
summary(O.ordata_OSOD_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.325993 0.674007
#mu     0.427758 2.260303
#sigma  0.283360 0.951113
#loglik at estimate:  -10651.58  

# Generate list of paralogs 

O.ordata_OSOD_WGD_paralogs <- O.ordata_OSOD  %>% 
  filter(Ks > 0.144398) %>% 
  filter(Ks < 0.711118) %>% 
  select(X)

write.table(O.ordata_OSOD_WGD_paralogs, quote = F, file = "Osmolindsaea_odorata_OSOD_WGD_paralogs.tsv")

#### Osmolindsaea ordata OSOR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.ordata_OSOR <- read.delim("ks_distributions/Osmolindsaea_odorata_OSOR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.ordata_OSOR_filt <- O.ordata_OSOR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.ordata_OSOR_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.ordata_OSOR_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Osmolindsaea ordata (OSOR)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.ordata_OSOR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.ordata_OSOR_filt.0 <- O.ordata_OSOR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.ordata_OSOR_normalmixEM <- normalmixEM(O.ordata_OSOR_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2 
summary(O.ordata_OSOR_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.420655 0.579345
#mu     0.425829 2.136353
#sigma  0.280219 0.980761
#loglik at estimate:  -10175.63 

# Generate list of paralogs 

O.ordata_OSOR_WGD_paralogs <- O.ordata_OSOR  %>% 
  filter(Ks > 0.14561) %>% 
  filter(Ks < 0.706048) %>% 
  select(X)

write.table(O.ordata_OSOR_WGD_paralogs, quote = F, file = "Osmolindsaea_odorata_OSOR_WGD_paralogs.tsv")

#### Omsunda japonia OSJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.japonica_OSJA <- read.delim("ks_distributions/Osmunda_japonica.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.japonica_OSJA_filt <- O.japonica_OSJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.japonica_OSJA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.japonica_OSJA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Osmunda japonica (OSJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.japonica_OSJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.japonica_OSJA_filt.0 <- O.japonica_OSJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.japonica_OSJA_normalmixEM <- normalmixEM(O.japonica_OSJA_filt.0$Ks, k=3)
# Fits to 1 WGD peak, ignore comp 2 
summary(O.japonica_OSJA_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.136620 0.383495 0.479885
#mu     0.160276 0.940986 2.605559
#sigma  0.113264 0.415542 0.756155
#loglik at estimate:  -8061.231 

# Generate list of paralogs 

O.japonica_OSJA_WGD_paralogs <- O.japonica_OSJA  %>% 
  filter(Ks > 0.525444) %>% 
  filter(Ks < 1.356528) %>% 
  select(X)

write.table(O.japonica_OSJA_WGD_paralogs, quote = F, file = "Osmunda_japonica_OSJA_WGD_paralogs.tsv")

#### Osmunda javanica VIBO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.javanica_VIBO <- read.delim("ks_distributions/Osmunda_javanica_VIBO.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.javanica_VIBO_filt <- O.javanica_VIBO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.javanica_VIBO_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.javanica_VIBO_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Osmunda javanica (VIBO)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.javanica_VIBO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.javanica_VIBO_filt.0 <- O.javanica_VIBO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.javanica_VIBO_normalmixEM <- normalmixEM(O.javanica_VIBO_filt.0$Ks, k=3)
# Fits to 1 WGD peak, ignore last component 
summary(O.javanica_VIBO_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.160652 0.402683 0.436665
#mu     0.128962 0.969304 2.622775
#sigma  0.106916 0.448995 0.775804
#loglik at estimate:  -6511.159

# Generate list of paralogs 

O.javanica_VIBO_WGD_paralogs <- O.javanica_VIBO  %>% 
  filter(Ks > 0.520309) %>% 
  filter(Ks < 1.418299) %>% 
  select(X)

write.table(O.javanica_VIBO_WGD_paralogs, quote = F, file = "Osmunda_javanica_VIBO_WGD_paralogs.tsv")

#### Osmunda sp. UOMY  ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

O.sp_UOMY <- read.delim("ks_distributions/Osmunda_sp._UOMY.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

O.sp_UOMY_filt <- O.sp_UOMY %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(O.sp_UOMY_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(O.sp_UOMY_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Osmunda sp. (UOMY)") + theme(plot.title = element_text(face = "italic"))

ggsave("O.sp_UOMY.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

O.sp_UOMY_filt.0 <- O.sp_UOMY_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

O.sp_UOMY_normalmixEM <- normalmixEM(O.sp_UOMY_filt.0$Ks, k=3)
# Fits to 1 WGD peak, ignore last component 
summary(O.sp_UOMY_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.174404 0.406717 0.418880
#mu     0.152871 1.012967 2.677388
#sigma  0.120852 0.450291 0.720565
#loglik at estimate:  -7207.229 

# Generate list of paralogs 

O.sp_UOMY_WGD_paralogs <- O.sp_UOMY  %>% 
  filter(Ks > 0.562676) %>% 
  filter(Ks < 1.463258) %>% 
  select(X)

write.table(O.sp_UOMY_WGD_paralogs, quote = F, file = "Osmunda_sp_UOMY_WGD_paralogs.tsv")

#### Phegopteris decursive-pinnata PHDP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.decursive_PHDP <- read.delim("ks_distributions/Phegopteris_ducursive-pinnata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.decursive_PHDP_filt <- P.decursive_PHDP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.decursive_PHDP_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.decursive_PHDP_filt, mapping = aes(x=Ks, ..scaled..*1300), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Phegopteris decursive-pinnata (PHDP)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.decursive_PHDP.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Phlebodium pseudoaureum ZQYU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.pseudoaureum_ZQYU <- read.delim("ks_distributions/Phlebodium_pseudoaureum_ZQYU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.pseudoaureum_ZQYU_filt <- P.pseudoaureum_ZQYU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.pseudoaureum_ZQYU_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.pseudoaureum_ZQYU_filt, mapping = aes(x=Ks, ..scaled..*1750), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Phlebodium pseudoaureum (ZQYU)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.pseudoaureum_ZQYU.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD


#### Phymatosorus grossus ORJE #### 

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.grossus_ORJE <- read.delim("ks_distributions/Phymatosorus_grossus_ORJE.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.grossus_ORJE_filt <- P.grossus_ORJE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.grossus_ORJE_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.grossus_ORJE_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Phymatosorus grossus (ORJE)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.grossus_ORJE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.grossus_ORJE_filt.0 <- P.grossus_ORJE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.grossus_ORJE_normalmixEM <- normalmixEM(P.grossus_ORJE_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore last component 
summary(P.grossus_ORJE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1972861 0.802714
#mu     0.0380227 1.950668
#sigma  0.0331745 1.081625
#loglik at estimate:  -6582.711

# Generate list of paralogs 

P.grossus_ORJE_WGD_paralogs <- P.grossus_ORJE  %>% 
  filter(Ks > 0.869043) %>% 
  filter(Ks < 3.032293) %>% 
  select(X)

write.table(P.grossus_ORJE_WGD_paralogs, quote = F, file = "Phymatosorus_grossus_ORJE_WGD_paralogs.tsv")

#### Pilularia globulifera KIIX ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.globulifera_KIIX <- read.delim("ks_distributions/Pilularia_globulifera_KIIX.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.globulifera_KIIX_filt <- P.globulifera_KIIX %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.globulifera_KIIX_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.globulifera_KIIX_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pilularia globulifera (KIIX)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.globulifera_KIIX.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.globulifera_KIIX_filt.0 <- P.globulifera_KIIX_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.globulifera_KIIX_normalmixEM <- normalmixEM(P.globulifera_KIIX_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.globulifera_KIIX_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1972861 0.802714
#mu     0.0380227 1.950668
#sigma  0.0331745 1.081625
#loglik at estimate:  -6582.711 

# Generate list of paralogs 

P.globulifera_KIIX_WGD_paralogs <- P.globulifera_KIIX  %>% 
  filter(Ks > 1.1224) %>% 
  filter(Ks < 3.223922) %>% 
  select(X)

write.table(P.globulifera_KIIX_WGD_paralogs, quote = F, file = "Pilularia_globulifera_KII_WGD_paralogs.tsv")
 
#### Pityrogramma trifoliata UJTT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.trifoliata_UJTT <- read.delim("ks_distributions/Pityrogramma_trifoliata_UJTT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.trifoliata_UJTT_filt <- P.trifoliata_UJTT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.trifoliata_UJTT_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.trifoliata_UJTT_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pityrogramma trifoliata (UJTT)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.trifoliata_UJTT.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Plagiogyria japonica UWOD #### 

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.japonica_UWOD <- read.delim("ks_distributions/Plagiogyria_japonica_UWOD.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.japonica_UWOD_filt <- P.japonica_UWOD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.japonica_UWOD_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.japonica_UWOD_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Plagiogyria japonica (UWOD)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.japonica_UWOD.png", height = 5, width = 8, dpi = 300)
 
# No evidence for WGD? 

#### Plagiogyria japonica PLJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.japonica_PLJA <- read.delim("ks_distributions/Plagiogyria_japonica_PLJA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.japonica_PLJA_filt <- P.japonica_PLJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.japonica_PLJA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.japonica_PLJA_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Plagiogyria japonica (PLJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.japonica_PLJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.japonica_PLJA_filt.0 <- P.japonica_PLJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.japonica_PLJA_normalmixEM <- normalmixEM(P.japonica_PLJA_filt.0$Ks, k=2)
# Fits to 2 WGD peaks 
summary(P.japonica_PLJA_normalmixEM)
#summary of normalmixEM object:
#  comp 1  comp 2
#lambda 0.323850 0.67615
#mu     0.206736 1.83120
#sigma  0.113262 1.00643
#loglik at estimate:  -8344.219

# Generate list of paralogs 

P.japonica_PLJA_WGD_paralogs_peak1 <- P.japonica_PLJA %>% 
  filter(Ks > 0.093474) %>% 
  filter(Ks < 0.319998) %>% 
  select(X)

write.table(P.japonica_PLJA_WGD_paralogs_peak1, quote = F, file = "P.japonica_PLJA_WGD_paralogs_peak1.tsv")

P.japonica_PLJA_WGD_paralogs_peak2 <- P.japonica_PLJA %>% 
  filter(Ks > 0.82477) %>% 
  filter(Ks < 2.8763) %>% 
  select(X)

write.table(P.japonica_PLJA_WGD_paralogs_peak2, quote = F, file = "P.japonica_PLJA_WGD_paralogs_peak2.tsv")

#### Plagiogyria stenoptera PLST ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.stenoptera_PLST <- read.delim("ks_distributions/Plagiogyria_stenoptera.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.stenoptera_PLST_filt <- P.stenoptera_PLST %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.stenoptera_PLST_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.stenoptera_PLST_filt, mapping = aes(x=Ks, ..scaled..*1400), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Plagiogyria stenoptera (PLST)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.stenoptera_PLST.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.stenoptera_PLST_filt.0 <- P.stenoptera_PLST_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.stenoptera_PLST_normalmixEM <- normalmixEM(P.stenoptera_PLST_filt.0$Ks, k=2)
# Fits to 2 WGD peaks? 
summary(P.stenoptera_PLST_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.244325 0.755675
#mu     0.202564 1.999983
#sigma  0.117396 1.010476
#loglik at estimate:  -11249.86 

# Generate list of paralogs 

P.stenoptera_PLST_WGD_paralogs_peak1 <- P.stenoptera_PLST %>% 
  filter(Ks > 0.085168) %>% 
  filter(Ks < 0.31996) %>% 
  select(X)

write.table(P.stenoptera_PLST_WGD_paralogs_peak1, quote = F, file = "P.stenoptera_PLST_WGD_paralogs_peak1.tsv")

P.stenoptera_PLST_WGD_paralogs_peak2 <- P.stenoptera_PLST %>% 
  filter(Ks > 0.989507) %>% 
  filter(Ks < 3.010459) %>% 
  select(X)

write.table(P.stenoptera_PLST_WGD_paralogs_peak2, quote = F, file = "P.stenoptera_PLST_WGD_paralogs_peak2.tsv")

#### Platycerium bifurcatum PLBI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.bifurcatum_PLBI <- read.delim("ks_distributions/Platycerium_bifurcatum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.bifurcatum_PLBI_filt <- P.bifurcatum_PLBI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.bifurcatum_PLBI_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.bifurcatum_PLBI_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Platycerium bifurcatum (PLBI)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.bifurcatum_PLBI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.bifurcatum_PLBI_filt.0 <- P.bifurcatum_PLBI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.bifurcatum_PLBI_normalmixEM <- normalmixEM(P.bifurcatum_PLBI_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.bifurcatum_PLBI_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.0955372 0.904463
#mu     0.1379522 2.218925
#sigma  0.1083563 0.979318
#loglik at estimate:  -4019.591

# Generate list of paralogs 

P.bifurcatum_PLBI_WGD_paralogs <- P.bifurcatum_PLBI %>% 
  filter(Ks > 1.239607) %>% 
  filter(Ks < 3.198243) %>% 
  select(X)

write.table(P.bifurcatum_PLBI_WGD_paralogs, quote = F, file = "Platycerium_bifurcatum_PLBI_WGD_paralogs.tsv")

#### Platycerium elephantotis PLEL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.elephantotis_PLEL <- read.delim("ks_distributions/Platycerium_elephantotis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.elephantotis_PLEL_filt <- P.elephantotis_PLEL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.elephantotis_PLEL_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.elephantotis_PLEL_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Platycerium elephantotis (PLEL)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.elephantotis_PLEL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.elephantotis_PLEL_filt.0 <- P.elephantotis_PLEL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.elephantotis_PLEL_normalmixEM <- normalmixEM(P.elephantotis_PLEL_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(P.elephantotis_PLEL_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.154909 0.845091
#mu     0.306253 2.351236
#sigma  0.223961 0.909021
#loglik at estimate:  -7326.091  

# Generate list of paralogs 

P.elephantotis_PLEL_WGD_paralogs <- P.elephantotis_PLEL %>% 
  filter(Ks > 1.442215) %>% 
  filter(Ks < 3.260257) %>% 
  select(X)

write.table(P.elephantotis_PLEL_WGD_paralogs, quote = F, file = "Platycerium_elephantotis_PLEL_WGD_paralogs.tsv")

#### Plenasium benksiifolium PLBA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.benksiifolium_PLBA <- read.delim("ks_distributions/Plenasium_banksiifolium.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.benksiifolium_PLBA_filt <- P.benksiifolium_PLBA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.benksiifolium_PLBA_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.benksiifolium_PLBA_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Plenasium benksiifolium (PLBA)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.benksiifolium_PLBA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.benksiifolium_PLBA_filt.0 <- P.benksiifolium_PLBA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.benksiifolium_PLBA_normalmixEM <- normalmixEM(P.benksiifolium_PLBA_filt.0$Ks, k=2)
# Fits to 1 WGD peak, comp 2 overfit? 
summary(P.benksiifolium_PLBA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.566648 0.433352
#mu     0.890700 2.750830
#sigma  0.575909 0.713614
#loglik at estimate:  -8333.354  

# Generate list of paralogs 

P.benksiifolium_PLBA_WGD_paralogs <- P.benksiifolium_PLBA %>% 
  filter(Ks > 0.314791) %>% 
  filter(Ks < 1.466609) %>% 
  select(X)

write.table(P.benksiifolium_PLBA_WGD_paralogs, quote = F, file = "Plenasium_benksiifolium_PLBA_WGD_paralogs.tsv")

#### Pleocnemia wintii PLWI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.wintii_PLWI <- read.delim("ks_distributions/Pleocnemia_winitii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.wintii_PLWI_filt <- P.wintii_PLWI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.wintii_PLWI_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.wintii_PLWI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pleocnemia wintii (PLWI)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.wintii_PLWI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.wintii_PLWI_filt.0 <- P.wintii_PLWI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.wintii_PLWI_normalmixEM <- normalmixEM(P.wintii_PLWI_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(P.wintii_PLWI_normalmixEM)
#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.128083 0.871917
#mu     0.208340 2.212477
#sigma  0.167517 0.965483
#loglik at estimate:  -8058.701   

# Generate list of paralogs 

P.wintii_PLWI_WGD_paralogs <- P.wintii_PLWI %>% 
  filter(Ks > 1.246994) %>% 
  filter(Ks < 3.17796) %>% 
  select(X)

write.table(P.wintii_PLWI_WGD_paralogs, quote = F, file = "Pleocnemia_winitii_PLWI_WGD_paralogs.tsv")

#### Pleopeltis polypodioides UJWU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.polypodioides_UJWU <- read.delim("ks_distributions/Pleopeltis_polypodioides_UJWU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.polypodioides_UJWU_filt <- P.polypodioides_UJWU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.polypodioides_UJWU_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.polypodioides_UJWU_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pleopeltis polypodioides (UJWU)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.polypodioides_UJWU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.polypodioides_UJWU_filt.0 <- P.polypodioides_UJWU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.polypodioides_UJWU_normalmixEM <- normalmixEM(P.polypodioides_UJWU_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.polypodioides_UJWU_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1656384 0.834362
#mu     0.1161239 2.190642
#sigma  0.0996032 1.003499
#loglik at estimate:  -4498.743   

# Generate list of paralogs 

P.polypodioides_UJWU_WGD_paralogs <- P.polypodioides_UJWU %>% 
  filter(Ks > 1.187143) %>% 
  filter(Ks < 3.194141) %>% 
  select(X)

write.table(P.polypodioides_UJWU_WGD_paralogs, quote = F, file = "Pleopeltis_polypodioides_UJWU_WGD_paralogs.tsv")

#### Pleurosoriopsis makinoi PLMA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.makinoi_PLMA <- read.delim("ks_distributions/Pleurosoriopsis_makinoi.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.makinoi_PLMA_filt <- P.makinoi_PLMA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.makinoi_PLMA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.makinoi_PLMA_filt, mapping = aes(x=Ks, ..scaled..*400), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pleurosoriopsis makinoi (PLMA)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.makinoi_PLMA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.makinoi_PLMA_filt.0 <- P.makinoi_PLMA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.makinoi_PLMA_normalmixEM <- normalmixEM(P.makinoi_PLMA_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.makinoi_PLMA_normalmixEM)
#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.140075 0.859925
#mu     0.292031 2.313279
#sigma  0.210092 0.932500
#loglik at estimate:  -5569.744 

# Generate list of paralogs 

P.makinoi_PLMA_WGD_paralogs <- P.makinoi_PLMA %>% 
  filter(Ks > 1.380779) %>% 
  filter(Ks < 3.245779) %>% 
  select(X)

write.table(P.makinoi_PLMA_WGD_paralogs, quote = F, file = "Pleurosoriopsis_makinoi_PLMA_WGD_paralogs.tsv")

#### Polypodium amorphum YLJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.amorphum_YLJA <- read.delim("ks_distributions/Polypodium_amorphum_YLJA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.amorphum_YLJA_filt <- P.amorphum_YLJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.amorphum_YLJA_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.amorphum_YLJA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polypodium amorphum (YLJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.amorphum_YLJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.amorphum_YLJA_filt.0 <- P.amorphum_YLJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.amorphum_YLJA_normalmixEM <- normalmixEM(P.amorphum_YLJA_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.amorphum_YLJA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.167751 0.832249
#mu     0.188513 2.105762
#sigma  0.155777 0.982773
#loglik at estimate:  -8117.524 

# Generate list of paralogs 

P.amorphum_YLJA_WGD_paralogs <- P.amorphum_YLJA %>% 
  filter(Ks > 1.122989) %>% 
  filter(Ks < 3.088535) %>% 
  select(X)

write.table(P.amorphum_YLJA_WGD_paralogs, quote = F, file = "Polypodium_amorphum_YLJA_WGD_paralogs.tsv")

#### Polypodium glycyrrhiza CJNT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.glycyrrhiza_CJNT <- read.delim("ks_distributions/Polypodium_glycyrrhiza_CJNT.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.glycyrrhiza_CJNT_filt <- P.glycyrrhiza_CJNT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.glycyrrhiza_CJNT_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.glycyrrhiza_CJNT_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polypodium glycyrrhiza (CJNT)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.glycyrrhiza_CJNT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.glycyrrhiza_CJNT_filt.0 <- P.glycyrrhiza_CJNT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.glycyrrhiza_CJNT_normalmixEM <- normalmixEM(P.glycyrrhiza_CJNT_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.glycyrrhiza_CJNT_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.158680 0.841320
#mu     0.119633 2.184610
#sigma  0.104802 0.992904
#loglik at estimate:  -7691.935 

# Generate list of paralogs 

P.glycyrrhiza_CJNT_WGD_paralogs <- P.glycyrrhiza_CJNT %>% 
  filter(Ks > 1.191706) %>% 
  filter(Ks < 3.177514) %>% 
  select(X)

write.table(P.glycyrrhiza_CJNT_WGD_paralogs, quote = F, file = "Polypodium_glycyrrhiza_CJNT_WGD_paralogs.tsv")

#### Polypodium hesperium GYFU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.hesperium_GYFU <- read.delim("ks_distributions/Polypodium_hesperium_GYFU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.hesperium_GYFU_filt <- P.hesperium_GYFU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.hesperium_GYFU_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.hesperium_GYFU_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polypodium hesperium (GYFU)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.hesperium_GYFU.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Polypodium hesperium IXLH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.hesperium_IXLH <- read.delim("ks_distributions/Polypodium_hesperium_IXLH.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.hesperium_IXLH_filt <- P.hesperium_IXLH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.hesperium_IXLH_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.hesperium_IXLH_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polypodium hesperium (IXLH)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.hesperium_IXLH.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Polypodium virginianum POVI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.virginianum_POVI <- read.delim("ks_distributions/Polypodium_virginianum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.virginianum_POVI_filt <- P.virginianum_POVI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.virginianum_POVI_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.virginianum_POVI_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polypodium virginianum (POVI)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.virginianum_POVI.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD?

#### Polystichum acrostichoides FQGQ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.acrostichoides_FQGQ <- read.delim("ks_distributions/Polystichum_acrostichoides_FQGQ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.acrostichoides_FQGQ_filt <- P.acrostichoides_FQGQ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.acrostichoides_FQGQ_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.acrostichoides_FQGQ_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polystichum acrostichoides (FQGQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.acrostichoides_FQGQ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.acrostichoides_FQGQ_filt.0 <- P.acrostichoides_FQGQ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.acrostichoides_FQGQ_normalmixEM <- normalmixEM(P.acrostichoides_FQGQ_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.acrostichoides_FQGQ_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.1064565 0.893543
#mu     0.1076636 2.076442
#sigma  0.0946859 1.004265
#loglik at estimate:  -6861.717 

# Generate list of paralogs 

P.acrostichoides_FQGQ_WGD_paralogs <- P.acrostichoides_FQGQ %>% 
  filter(Ks > 1.072177) %>% 
  filter(Ks < 3.080707) %>% 
  select(X)

write.table(P.acrostichoides_FQGQ_WGD_paralogs, quote = F, file = "Polystichum_acrostichoides_FQGQ_WGD_paralogs.tsv")

#### Polystichum tripteron POTR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.tripteron_POTR <- read.delim("ks_distributions/Polystichum_tripteron.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.tripteron_POTR_filt <- P.tripteron_POTR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.tripteron_POTR_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.tripteron_POTR_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Polystichum tripteron (POTR)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.tripteron_POTR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.tripteron_POTR_filt.0 <- P.tripteron_POTR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.tripteron_POTR_normalmixEM <- normalmixEM(P.tripteron_POTR_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(P.tripteron_POTR_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.108428 0.891572
#mu     0.190059 2.144853
#sigma  0.145302 0.963123
#loglik at estimate:  -7694.074  

# Generate list of paralogs 

P.tripteron_POTR_WGD_paralogs <- P.tripteron_POTR %>% 
  filter(Ks > 1.18173) %>% 
  filter(Ks < 3.107976) %>% 
  select(X)

write.table(P.tripteron_POTR_WGD_paralogs, quote = F, file = "Polystichum_tripteron_POTR_WGD_paralogs.tsv")

#### Pronephrium simplex PRSI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.simplex_PRSI <- read.delim("ks_distributions/Pronephrium_simplex.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.simplex_PRSI_filt <- P.simplex_PRSI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.simplex_PRSI_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.simplex_PRSI_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pronephrium simplex (PRSI)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.simplex_PRSI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.simplex_PRSI_filt.0 <- P.simplex_PRSI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.simplex_PRSI_normalmixEM <- normalmixEM(P.simplex_PRSI_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 1
summary(P.simplex_PRSI_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.184692 0.815308
#mu     3.423823 1.630355
#sigma  0.348364 0.912791
#loglik at estimate:  -5342.636 

# Generate list of paralogs 

P.simplex_PRSI_WGD_paralogs <- P.simplex_PRSI %>% 
  filter(Ks > 0.717564) %>% 
  filter(Ks < 2.543146) %>% 
  select(X)

write.table(P.simplex_PRSI_WGD_paralogs, quote = F, file = "Pronephrium_simplex_PRSI_WGD_paralogs.tsv")

#### Prosaptia obliquata PROB ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.obliquata_PROB <- read.delim("ks_distributions/Prosaptia_obliquata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.obliquata_PROB_filt <- P.obliquata_PROB %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.obliquata_PROB_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.obliquata_PROB_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Prosaptia obliquata (PROB)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.obliquata_PROB.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Psilotum nudum QVMR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.nudum_QVMR <- read.delim("ks_distributions/Psilotum_nudum_QVMR.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.nudum_QVMR_filt <- P.nudum_QVMR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.nudum_QVMR_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.nudum_QVMR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Psilotum nudum (QVMR)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.nudum_QVMR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.nudum_QVMR_filt.0 <- P.nudum_QVMR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.nudum_QVMR_normalmixEM <- normalmixEM(P.nudum_QVMR_filt.0$Ks, k=5,maxit = 2000)
# Fits to 2 WGD peak- comps 2 and 3 
summary(P.nudum_QVMR_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3   comp 4   comp 5
#lambda 0.0806818 0.271179 0.304280 0.294062 0.049797
#mu     0.0360708 0.376602 1.063551 2.438207 3.629495
#sigma  0.0299379 0.172657 0.403566 0.634578 0.215744
#loglik at estimate:  -3949.928

# Generate list of paralogs 

P.nudum_QVMR_WGD_paralogs_peak1 <- P.nudum_QVMR %>% 
  filter(Ks > 0.143479) %>% 
  filter(Ks < 0.539893) %>% 
  select(X)

write.table(P.nudum_QVMR_WGD_paralogs_peak1, quote = F, file = "Psilotum_nudum_QVMR_WGD_paralogs_peak1.tsv")

P.nudum_QVMR_WGD_paralogs_peak2 <- P.nudum_QVMR %>% 
  filter(Ks > 0.683146) %>% 
  filter(Ks < 1.717604) %>% 
  select(X)

write.table(P.nudum_QVMR_WGD_paralogs_peak2, quote = F, file = "Psilotum_nudum_QVMR_WGD_paralogs_peak2.tsv")

#### Psilotum nudum PSNU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.nudum_PSNU <- read.delim("ks_distributions/Psilotum_nudum_PSNU.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.nudum_PSNU_filt <- P.nudum_PSNU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.nudum_PSNU_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.nudum_PSNU_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Psioltum nudum (PSNU)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.nudum_PSNU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.nudum_PSNU_filt.0 <- P.nudum_PSNU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.nudum_PSNU_normalmixEM <- normalmixEM(P.nudum_PSNU_filt.0$Ks, k=3, maxit = 2000)
# Fits to 2 WGD peaks, ignore comp 3
summary(P.nudum_PSNU_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.287293 0.405024 0.307683
#mu     0.255840 1.110968 2.802823
#sigma  0.196937 0.503737 0.674483
#loglik at estimate:  -8140.5 

# Generate list of paralogs 

P.nudum_PSNU_WGD_paralogs_peak1 <- P.nudum_PSNU %>% 
  filter(Ks > 0.058903) %>% 
  filter(Ks < 0.452777) %>% 
  select(X)

write.table(P.nudum_PSNU_WGD_paralogs_peak1, quote = F, file = "Psilotum_nudum_PSNU_WGD_paralogs_peak1.tsv")

P.nudum_PSNU_WGD_paralogs_peak2 <- P.nudum_PSNU %>% 
  filter(Ks > 0.607231) %>% 
  filter(Ks < 1.614705) %>% 
  select(X)

write.table(P.nudum_PSNU_WGD_paralogs_peak2, quote = F, file = "Psilotum_nudum_PSNU_WGD_paralogs_peak2.tsv")

#### Psioltum nudum PSND ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.nudum_PSND <- read.delim("ks_distributions/Psilotum_nudum_PSND.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.nudum_PSND_filt <- P.nudum_PSND %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.nudum_PSND_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.nudum_PSND_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Psioltum nudum (PSND)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.nudum_PSND.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.nudum_PSND_filt.0 <- P.nudum_PSND_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.nudum_PSND_normalmixEM <- normalmixEM(P.nudum_PSND_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3
summary(P.nudum_PSND_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.311049 0.431045 0.257906
#mu     0.341686 1.200375 2.863526
#sigma  0.198207 0.517229 0.653415
#loglik at estimate:  -8515.928 

# Generate list of paralogs 

P.nudum_PSND_WGD_paralogs_peak1 <- P.nudum_PSND %>% 
  filter(Ks > 0.203945) %>% 
  filter(Ks < 0.549259) %>% 
  select(X)

write.table(P.nudum_PSND_WGD_paralogs_peak1, quote = F, file = "Psilotum_nudum_PSND_WGD_paralogs_peak1.tsv")

P.nudum_PSND_WGD_paralogs_peak2 <- P.nudum_PSND %>% 
  filter(Ks > 0.659985) %>% 
  filter(Ks < 1.467117) %>% 
  select(X)

write.table(P.nudum_PSND_WGD_paralogs_peak2, quote = F, file = "Psilotum_nudum_PSND_WGD_paralogs_peak2.tsv")

#### Pteridium aquilinum subsp. latiusculum PTAQ #### 

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.aquilinum_PTAQ <- read.delim("ks_distributions/Pteridium_aquilinum_subsp._latisculum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.aquilinum_PTAQ_filt <- P.aquilinum_PTAQ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.aquilinum_PTAQ_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.aquilinum_PTAQ_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteridium aquilinum subsp. latiusculum (PTAQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.aquilinum_PTAQ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.aquilinum_PTAQ_filt.0 <- P.aquilinum_PTAQ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.aquilinum_PTAQ_normalmixEM <- normalmixEM(P.aquilinum_PTAQ_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(P.aquilinum_PTAQ_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.752445 0.247555
#mu     1.469501 3.214400
#sigma  0.871570 0.477275
#loglik at estimate:  -9330.993 

# Generate list of paralogs 

P.aquilinum_PTAQ_WGD_paralogs <- P.aquilinum_PTAQ %>% 
  filter(Ks > 0.597931) %>% 
  filter(Ks < 2.341071) %>% 
  select(X)

write.table(P.aquilinum_PTAQ_WGD_paralogs, quote = F, file = "Pteridium_aquilinum_PTAQ_WGD_paralogs.tsv")

#### Pteridium revolutum PTRE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.revolutum_PTRE <- read.delim("ks_distributions/Pteridium_revolutum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.revolutum_PTRE_filt <- P.revolutum_PTRE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.revolutum_PTRE_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.revolutum_PTRE_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteridium revolutum (PTRE)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.revolutum_PTRE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.revolutum_PTRE_filt.0 <- P.revolutum_PTRE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.revolutum_PTRE_normalmixEM <- normalmixEM(P.revolutum_PTRE_filt.0$Ks, k=2, maxit = 2000)
# Fits to 2 WGD peaks, ignore comp 3
summary(P.revolutum_PTRE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0914185 0.908582
#mu     0.1718078 2.076863
#sigma  0.1139685 0.979440
#loglik at estimate:  -5578.469 

# Generate list of paralogs 

P.revolutum_PTRE_WGD_paralogs <- P.revolutum_PTRE %>% 
  filter(Ks > 1.097423) %>% 
  filter(Ks < 3.056303) %>% 
  select(X)

write.table(P.revolutum_PTRE_WGD_paralogs, quote = F, file = "Pteridium_revolutum_PTRE_WGD_paralogs.tsv")

#### Pteridrys cnemidaria PTCN ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.cnemidaria_PTCN <- read.delim("ks_distributions/Pteridrys_cnemidaria.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.cnemidaria_PTCN_filt <- P.cnemidaria_PTCN %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.cnemidaria_PTCN_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.cnemidaria_PTCN_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteridrys cnemidaria (PTCN)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.cnemidaria_PTCN.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.cnemidaria_PTCN_filt.0 <- P.cnemidaria_PTCN_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.cnemidaria_PTCN_normalmixEM <- normalmixEM(P.cnemidaria_PTCN_filt.0$Ks, k=2, maxit = 2000)
# Fits to 2 WGD peaks, ignore comp 3
summary(P.cnemidaria_PTCN_normalmixEM)
#  summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.106971 0.893029
#mu     0.139950 2.133222
#sigma  0.121391 0.999808
#loglik at estimate:  -9009.984  

# Generate list of paralogs 

P.cnemidaria_PTCN_WGD_paralogs <- P.cnemidaria_PTCN %>% 
  filter(Ks > 1.133414) %>% 
  filter(Ks < 3.13303) %>% 
  select(X)

write.table(P.cnemidaria_PTCN_WGD_paralogs, quote = F, file = "Pteridrys_cnemidaria_PTCN_WGD_paralogs.tsv")

#### Pteris ensiformis FLTD ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.ensiformis_FLTD <- read.delim("ks_distributions/Pteris_ensiformis_FLTD.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.ensiformis_FLTD_filt <- P.ensiformis_FLTD %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.ensiformis_FLTD_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.ensiformis_FLTD_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteris ensiformis (FLTD)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.ensiformis_FLTD.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.ensiformis_FLTD_filt.0 <- P.ensiformis_FLTD_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.ensiformis_FLTD_normalmixEM <- normalmixEM(P.ensiformis_FLTD_filt.0$Ks, k=2)
# Fits to 1 WGD
summary(P.ensiformis_FLTD_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1759885 0.824011
#mu     0.0607004 2.088269
#sigma  0.0587028 1.050039
#loglik at estimate:  -4808.314 

# Generate list of paralogs 

P.ensiformis_FLTD_WGD_paralogs <- P.ensiformis_FLTD %>% 
  filter(Ks > 1.03823) %>% 
  filter(Ks < 3.138308) %>% 
  select(X)

write.table(P.ensiformis_FLTD_WGD_paralogs, quote = F, file = "Pteris_ensiformis_FLTD_WGD_paralogs.tsv")

#### Pteris fauriei PTFA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.fauriei_PTFA <- read.delim("ks_distributions/Pteris_fauriei.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.fauriei_PTFA_filt <- P.fauriei_PTFA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.fauriei_PTFA_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.fauriei_PTFA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteris fauriei (PTFA)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.fauriei_PTFA.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Pteris vittata POPJ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.vittata_POPJ <- read.delim("ks_distributions/Pteris_vittata_POPJ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.vittata_POPJ_filt <- P.vittata_POPJ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.vittata_POPJ_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.vittata_POPJ_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteris vittata (POPJ)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.vittata_POPJ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.vittata_POPJ_filt.0 <- P.vittata_POPJ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.vittata_POPJ_normalmixEM <- normalmixEM(P.vittata_POPJ_filt.0$Ks, k=2)
# Fits to 1 WGD
summary(P.vittata_POPJ_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.1572429 0.842757
#mu     0.0495594 2.106530
#sigma  0.0486652 1.047755
#loglik at estimate:  -6975.494 

# Generate list of paralogs 

P.vittata_POPJ_WGD_paralogs <- P.vittata_POPJ %>% 
  filter(Ks > 1.058775) %>% 
  filter(Ks < 3.154285) %>% 
  select(X)

write.table(P.vittata_POPJ_WGD_paralogs, quote = F, file = "Pteris_vittata_POPJ_WGD_paralogs.tsv")

#### Pteris vittata PTVI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.vittata_PTVI <- read.delim("ks_distributions/Pteris_vittata_PTVI.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.vittata_PTVI_filt <- P.vittata_PTVI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.vittata_PTVI_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.vittata_PTVI_filt, mapping = aes(x=Ks, ..scaled..*2000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pteris vittata (PTVI)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.vittata_PTVI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.vittata_PTVI_filt.0 <- P.vittata_PTVI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.vittata_PTVI_normalmixEM <- normalmixEM(P.vittata_PTVI_filt.0$Ks, k=2)
# Fits to 2 WGD peaks, ignore comp 3
summary(P.vittata_PTVI_normalmixEM)
#summary of normalmixEM object:
#           comp 1   comp 2
#lambda 0.1514880 0.848512
#mu     0.0771992 2.262044
#sigma  0.0762233 0.981208
#loglik at estimate:  -6500.541 

# Generate list of paralogs 

P.vittata_PTVI_WGD_paralogs <- P.vittata_PTVI %>% 
  filter(Ks > 1.280836) %>% 
  filter(Ks < 3.243252) %>% 
  select(X)

write.table(P.vittata_PTVI_WGD_paralogs, quote = F, file = "Pteris_vittata_PTVI_WGD_paralogs.tsv")

#### Ptisana pellucida PTPE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.pellucida_PTPE <- read.delim("ks_distributions/Ptisana_pellucida.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.pellucida_PTPE_filt <- P.pellucida_PTPE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.pellucida_PTPE_filt.0, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.pellucida_PTPE_filt.0, mapping = aes(x=Ks, ..scaled..*3500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Ptisana pellucida (PTPE)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.pellucida_PTPE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.pellucida_PTPE_filt.0 <- P.pellucida_PTPE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.pellucida_PTPE_normalmixEM <- normalmixEM(P.pellucida_PTPE_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD, comp 2 overfit?  
summary(P.pellucida_PTPE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.537766 0.462234
#mu     0.783637 2.459798
#sigma  0.433982 0.826722
#loglik at estimate:  -8656.165 

# Generate list of paralogs 

P.pellucida_PTPE_WGD_paralogs <- P.pellucida_PTPE %>% 
  filter(Ks > 0.349655) %>% 
  filter(Ks < 1.217619) %>% 
  select(X)

write.table(P.pellucida_PTPE_WGD_paralogs, quote = F, file = "Ptisana_pellucida_PTPE_WGD_paralogs.tsv")

#### Pyrrosia subfurfuracea PYSU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.subfurfuracea_PYSU <- read.delim("ks_distributions/Pyrrosia_subfuracea.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.subfurfuracea_PYSU_filt <- P.subfurfuracea_PYSU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.subfurfuracea_PYSU_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.subfurfuracea_PYSU_filt, mapping = aes(x=Ks, ..scaled..*3000), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Pyrrosia subfurfuracea (PYSU)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.subfurfuracea_PYSU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.subfurfuracea_PYSU_filt.0 <- P.subfurfuracea_PYSU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.subfurfuracea_PYSU_normalmixEM <- normalmixEM(P.subfurfuracea_PYSU_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD, comp 2 overfit?  
summary(P.subfurfuracea_PYSU_normalmixEM)
#summary of normalmixEM object:
#       comp 1   comp 2
#lambda 0.533259 0.466741
#mu     0.253915 2.095774
#sigma  0.139154 0.970664
#loglik at estimate:  -8844.638

# Generate list of paralogs 

P.subfurfuracea_PYSU_WGD_paralogs <- P.subfurfuracea_PYSU %>% 
  filter(Ks > 0.114761) %>% 
  filter(Ks < 0.393069) %>% 
  select(X)

write.table(P.subfurfuracea_PYSU_WGD_paralogs, quote = F, file = "Pyrrosia_subfurfuracea_PYSU_WGD_paralogs.tsv")

#### Rhachidosorus sp. RHSP ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

R.sp_RHSP <- read.delim("ks_distributions/Rhachidosorus_sp._XQ-2018.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

R.sp_RHSP_filt <- R.sp_RHSP %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(R.sp_RHSP_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(R.sp_RHSP_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Rhachidosorus sp. (RHSP)") + theme(plot.title = element_text(face = "italic"))

ggsave("R.sp_RHSP.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

R.sp_RHSP_filt.0 <- R.sp_RHSP_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

R.sp_RHSP_normalmixEM <- normalmixEM(R.sp_RHSP_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(R.sp_RHSP_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0954421 0.904558
#mu     0.1152511 2.066358
#sigma  0.0934433 0.993946
#loglik at estimate:  -8111.469

# Generate list of paralogs 

R.sp_RHSP_WGD_paralogs <- R.sp_RHSP %>% 
  filter(Ks > 1.072412) %>% 
  filter(Ks < 3.060304) %>% 
  select(X)

write.table(R.sp_RHSP_WGD_paralogs, quote = F, file = "Rhachidosorous_sp_RHSP_WGD_paralogs.tsv")

#### Rhachidosorus mesosorus RHME ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

P.mesosorus_RHME <- read.delim("ks_distributions/Rhachidosorus_mesosorus.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

P.mesosorus_RHME_filt <- P.mesosorus_RHME %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(P.mesosorus_RHME_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(P.mesosorus_RHME_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Rhachidosorus mesosorus (RHME)") + theme(plot.title = element_text(face = "italic"))

ggsave("P.mesosorus_RHME.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

P.mesosorus_RHME_filt.0 <- P.mesosorus_RHME_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

P.mesosorus_RHME_normalmixEM <- normalmixEM(P.mesosorus_RHME_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(P.mesosorus_RHME_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.827941 0.172059
#mu     1.431501 3.367928
#sigma  0.882472 0.365942
#loglik at estimate:  -9568.251 

# Generate list of paralogs 

P.mesosorus_RHME_WGD_paralogs <- P.mesosorus_RHME %>% 
  filter(Ks > 0.549029) %>% 
  filter(Ks < 2.313973) %>% 
  select(X)

write.table(P.mesosorus_RHME_WGD_paralogs, quote = F, file = "Rhachidosorous_mesosorus_RHME_WGD_paralogs.tsv")

#### Rhachidororus pulcher RHPU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

R.pulcher_RHPU <- read.delim("ks_distributions/Rhadchidosorus_pulcher.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

R.pulcher_RHPU_filt <- R.pulcher_RHPU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(R.pulcher_RHPU_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(R.pulcher_RHPU_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Rhachidosorus pulcher (RHPU)") + theme(plot.title = element_text(face = "italic"))

ggsave("R.pulcher_RHPU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

R.pulcher_RHPU_filt.0 <- R.pulcher_RHPU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

R.pulcher_RHPU_normalmixEM <- normalmixEM(R.pulcher_RHPU_filt.0$Ks, k=2)
# Fits to 1 WGD peak
summary(R.pulcher_RHPU_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0913745 0.908626
#mu     0.1586833 2.107731
#sigma  0.1300797 0.978991
#loglik at estimate:  -7984.547 

# Generate list of paralogs 

R.pulcher_RHPU_WGD_paralogs <- R.pulcher_RHPU %>% 
  filter(Ks > 1.12874) %>% 
  filter(Ks < 3.086722) %>% 
  select(X)

write.table(R.pulcher_RHPU_WGD_paralogs, quote = F, file = "Rhachidosorous_pulcher_RHPU_WGD_paralogs.tsv")

#### Salvinia natans SANT ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.natans_SANT <- read.delim("ks_distributions/Salvinia_natans.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.natans_SANT_filt <- S.natans_SANT %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.natans_SANT_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.natans_SANT_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Salvinia natans (SANT)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.natans_SANT.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.natans_SANT_filt.0 <- S.natans_SANT_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.natans_SANT_normalmixEM <- normalmixEM(S.natans_SANT_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2
summary(S.natans_SANT_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.110840 0.889160
#mu     0.112753 2.218075
#sigma  0.102629 0.989913
#loglik at estimate:  -6163.373  

# Generate list of paralogs 

S.natans_SANT_WGD_paralogs <- S.natans_SANT %>% 
  filter(Ks > 1.228162) %>% 
  filter(Ks < 3.207988) %>% 
  select(X)

write.table(S.natans_SANT_WGD_paralogs, quote = F, file = "Salvinia_natans_SANT_WGD_paralogs.tsv")

#### Salvinia natans SANA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.natans_SANA <- read.delim("ks_distributions/Salvinia_natans.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.natans_SANA_filt <- S.natans_SANA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.natans_SANA_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.natans_SANA_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Salvinia natans (SANA)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.natans_SANA.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Sceptridium dissectum EEAQ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.dissectum_EEAQ <- read.delim("ks_distributions/Sceptridium_dissectum_EEAQ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.dissectum_EEAQ_filt <- S.dissectum_EEAQ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.dissectum_EEAQ_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.dissectum_EEAQ_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Sceptridium dissectum (EEAQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.dissectum_EEAQ.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.dissectum_EEAQ_filt.0 <- S.dissectum_EEAQ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.dissectum_EEAQ_normalmixEM <- normalmixEM(S.dissectum_EEAQ_filt.0$Ks, k=3, maxit = 2000)
# Fits to 2 WGDs, ignore third component 
summary(S.dissectum_EEAQ_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2   comp 3
#lambda 0.424027 0.321191 0.254781
#mu     0.246341 1.153272 2.796617
#sigma  0.150430 0.499949 0.658370
#loglik at estimate:  -7204.033

# Generate list of paralogs 

S.dissectum_EEAQ_WGD_paralogs_peak1 <- S.dissectum_EEAQ %>% 
  filter(Ks > 0.095911) %>% 
  filter(Ks < 0.396771) %>% 
  select(X)

write.table(S.dissectum_EEAQ_WGD_paralogs_peak1, quote = F, file = "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak1.tsv")

S.dissectum_EEAQ_WGD_paralogs_peak2 <- S.dissectum_EEAQ %>% 
  filter(Ks > 0.653323) %>% 
  filter(Ks < 1.653221) %>% 
  select(X)

write.table(S.dissectum_EEAQ_WGD_paralogs_peak2, quote = F, file = "Sceptridium_dissectum_EEAQ_WGD_paralogs_peak2.tsv")

#### Sceptridium japonicum SCJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.japonicum_SCJA <- read.delim("ks_distributions/Sceptridium_japonicum.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.japonicum_SCJA_filt <- S.japonicum_SCJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.japonicum_SCJA_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.japonicum_SCJA_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Sceptridium japonicum (SCJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.japonicum_SCJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.japonicum_SCJA_filt.0 <- S.japonicum_SCJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.japonicum_SCJA_normalmixEM <- normalmixEM(S.japonicum_SCJA_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3
summary(S.japonicum_SCJA_normalmixEM)
#summary of normalmixEM object:
#       comp 1   comp 2   comp 3
#lambda 0.372493 0.360454 0.267053
#mu     0.272213 1.222062 2.841623
#sigma  0.159049 0.516639 0.648434
#loglik at estimate:  -7829.517  

# Generate list of paralogs 

S.japonicum_SCJA_WGD_paralogs_peak1 <- S.japonicum_SCJA %>% 
  filter(Ks > 0.113164) %>% 
  filter(Ks < 0.431262) %>% 
  select(X)

write.table(S.japonicum_SCJA_WGD_paralogs_peak1, quote = F, file = "Sceptridium_japonicum_SCJA_WGD_paralogs_peak1.tsv")

S.japonicum_SCJA_WGD_paralogs_peak2 <- S.japonicum_SCJA %>% 
  filter(Ks > 0.705423) %>% 
  filter(Ks < 1.738701) %>% 
  select(X)

write.table(S.japonicum_SCJA_WGD_paralogs_peak2, quote = F, file = "Sceptridium_japonicum_SCJA_WGD_paralogs_peak2.tsv")

#### Schizaea dichotomsa SCDI ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.dichotoma_SCDI <- read.delim("ks_distributions/Schizea_dichotoma.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.dichotoma_SCDI_filt <- S.dichotoma_SCDI %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.dichotoma_SCDI_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.dichotoma_SCDI_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Schizaea dichotoma (SCDI)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.dichotoma_SCDI.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.dichotoma_SCDI_filt.0 <- S.dichotoma_SCDI_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.dichotoma_SCDI_normalmixEM <- normalmixEM(S.dichotoma_SCDI_filt.0$Ks, k=2)
# Fits to 1 WGD peak, ignore comp 2? 
summary(S.dichotoma_SCDI_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.390357 0.609643
#mu     0.499099 2.239146
#sigma  0.279555 0.937293
#loglik at estimate:  -9494.379 

# Generate list of paralogs 

S.dichotoma_SCDI_WGD_paralogs <- S.dichotoma_SCDI %>% 
  filter(Ks > 0.219544) %>% 
  filter(Ks < 0.778654) %>% 
  select(X)

write.table(S.dichotoma_SCDI_WGD_paralogs, quote = F, file = "Schizaea_dichotoma_SCDI_WGD_paralogs.tsv")

#### Selliguea feei SEFE ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.feei_SEFE <- read.delim("ks_distributions/Selliguea_feei.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.feei_SEFE_filt <- S.feei_SEFE %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.feei_SEFE_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.feei_SEFE_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Selliguea feei (SEFE)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.feei_SEFE.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.feei_SEFE_filt.0 <- S.feei_SEFE_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.feei_SEFE_normalmixEM <- normalmixEM(S.feei_SEFE_filt.0$Ks, k=2)
# Fits to 1 WGD peak 
summary(S.feei_SEFE_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.146362 0.853638
#mu     0.234946 2.324360
#sigma  0.178017 0.922140
#loglik at estimate:  -8168.679 

# Generate list of paralogs 

S.feei_SEFE_WGD_paralogs <- S.feei_SEFE %>% 
  filter(Ks > 1.40222) %>% 
  filter(Ks < 3.2465) %>% 
  select(X)

write.table(S.feei_SEFE_WGD_paralogs, quote = F, file = "Selliguea_feei_SEFE_WGD_paralogs.tsv")

#### Stenochlaena palustris STPL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.palustris_STPL <- read.delim("ks_distributions/Stenochlaena_palustris_STPL.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.palustris_STPL_filt <- S.palustris_STPL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.palustris_STPL_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.palustris_STPL_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Stenochlaena palustris (STPL)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.palustris_STPL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.palustris_STPL_filt.0 <- S.palustris_STPL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.palustris_STPL_normalmixEM <- normalmixEM(S.palustris_STPL_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3 
summary(S.palustris_STPL_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.417453 0.480899 0.101649
#mu     0.335652 1.964659 3.467961
#sigma  0.157333 0.774009 0.323112
#loglik at estimate:  -8838.566

# Generate list of paralogs 

S.palustris_STPL_WGD_paralogs_peak1 <- S.palustris_STPL %>% 
  filter(Ks > 0.178319) %>% 
  filter(Ks < 0.492985) %>% 
  select(X)

write.table(S.palustris_STPL_WGD_paralogs_peak1, quote = F, file = "Stecnochlaena_palustris_STPL_WGD_paralogs_peak1.tsv")

S.palustris_STPL_WGD_paralogs_peak2 <- S.palustris_STPL %>% 
  filter(Ks > 1.19065) %>% 
  filter(Ks < 2.738668) %>% 
  select(X)

write.table(S.palustris_STPL_WGD_paralogs_peak2, quote = F, file = "Stecnochlaena_palustris_STPL_WGD_paralogs_peak2.tsv")

#### Stenochlaena palustris STPA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.palustris_STPA <- read.delim("ks_distributions/Stenochlaena_palustris_STPA.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.palustris_STPA_filt <- S.palustris_STPA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.palustris_STPA_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.palustris_STPA_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Stenochlaena palustris (STPA)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.palustris_STPA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.palustris_STPA_filt.0 <- S.palustris_STPA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.palustris_STPA_normalmixEM <- normalmixEM(S.palustris_STPA_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, comp 3 is overfit? 
summary(S.palustris_STPA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2    comp 3
#lambda 0.418362 0.487477 0.0941603
#mu     0.326528 1.919845 3.4645373
#sigma  0.165129 0.797285 0.3255461
#loglik at estimate:  -8754.305 

# Generate list of paralogs 

S.palustris_STPA_WGD_paralogs_peak1 <- S.palustris_STPA %>% 
  filter(Ks > 0.161399) %>% 
  filter(Ks < 0.491657) %>% 
  select(X)

write.table(S.palustris_STPA_WGD_paralogs_peak1, quote = F, file = "Stecnochlaena_palustris_STPA_WGD_paralogs_peak1.tsv")

S.palustris_STPA_WGD_paralogs_peak2 <- S.palustris_STPA %>% 
  filter(Ks > 1.12256) %>% 
  filter(Ks < 2.71713) %>% 
  select(X)

write.table(S.palustris_STPA_WGD_paralogs_peak2, quote = F, file = "Stecnochlaena_palustris_STPA_WGD_paralogs_peak2.tsv")

#### Stenograma wilfordii STFO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.wilfordii_STFO <- read.delim("ks_distributions/Stenogramma_wilfordii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.wilfordii_STFO_filt <- S.wilfordii_STFO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.wilfordii_STFO_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.wilfordii_STFO_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Stenogramma wilfordii (STFO)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.wilfordii_STFO.png", height = 5, width = 8, dpi = 300)

#No evidence of WGD? 

#### Sticherus truncatus STTR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

S.truncatus_STTR <- read.delim("ks_distributions/Sticherus_truncatus.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

S.truncatus_STTR_filt <- S.truncatus_STTR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(S.truncatus_STTR_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(S.truncatus_STTR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Sticherus truncatus (STTR)") + theme(plot.title = element_text(face = "italic"))

ggsave("S.truncatus_STTR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

S.truncatus_STTR_filt.0 <- S.truncatus_STTR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

S.truncatus_STTR_normalmixEM <- normalmixEM(S.truncatus_STTR_filt.0$Ks, k=2)
# Fits to 1 WGD, comp 1 is overfit? 
summary(S.truncatus_STTR_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.469498 0.362237 0.168264
#mu     1.986187 0.388932 3.379622
#sigma  0.762160 0.211392 0.364422
#loglik at estimate:  -7453.816 

# Generate list of paralogs 

S.truncatus_STTR_WGD_paralogs <- S.truncatus_STTR %>% 
  filter(Ks > 0.299014) %>% 
  filter(Ks < 1.172188) %>% 
  select(X)

write.table(S.truncatus_STTR_WGD_paralogs, quote = F, file = "Sticherus_truncatus_STTR_WGD_paralogs_peak1.tsv")

#### Taenitis blechnoides TABL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.blechnoides_TABL <- read.delim("ks_distributions/Taenitis_blechnoides.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.blechnoides_TABL_filt <- T.blechnoides_TABL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.blechnoides_TABL_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.blechnoides_TABL_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Taenitis blechnoides (TABL)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.blechnoides_TABL.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.blechnoides_TABL_filt.0 <- T.blechnoides_TABL_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.blechnoides_TABL_normalmixEM <- normalmixEM(T.blechnoides_TABL_filt.0$Ks, k=3)
# Fits to 2 WGD peaks, ignore comp 3 
summary(T.blechnoides_TABL_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.469498 0.362237 0.168264
#mu     1.986187 0.388932 3.379622
#sigma  0.762160 0.211392 0.364422
#loglik at estimate:  -7453.816 

# Generate list of paralogs 

T.blechnoides_TABL_WGD_paralogs_peak1 <- T.blechnoides_TABL %>% 
  filter(Ks > 0.17754) %>% 
  filter(Ks < 0.600324) %>% 
  select(X)

write.table(T.blechnoides_TABL_WGD_paralogs_peak1, quote = F, file = "Taenitis_blechnoides_TABL_WGD_paralogs_peak1.tsv")

T.blechnoides_TABL_WGD_paralogs_peak2 <- T.blechnoides_TABL %>% 
  filter(Ks > 1.224027) %>% 
  filter(Ks < 2.748347) %>% 
  select(X)

write.table(T.blechnoides_TABL_WGD_paralogs_peak2, quote = F, file = "Taenitis_blechnoides_TABL_WGD_paralogs_peak2.tsv")

#### Tectaria moresi TEMO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.moresi_TEMO <- read.delim("ks_distributions/Tectaria_morsei.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.moresi_TEMO_filt <- T.moresi_TEMO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.moresi_TEMO_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.moresi_TEMO_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Tectaria moresi (TEMO)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.moresi_TEMO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.moresi_TEMO_filt.0 <- T.moresi_TEMO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.moresi_TEMO_normalmixEM <- normalmixEM(T.moresi_TEMO_filt.0$Ks, k=2)
# Fits to 2 WGD peaks, ignore comp 3 
summary(T.moresi_TEMO_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.124373 0.875627
#mu     0.167616 2.117900
#sigma  0.123767 0.977425
#loglik at estimate:  -9438.401 

# Generate list of paralogs 

T.moresi_TEMO_WGD_paralogs <- T.moresi_TEMO %>% 
  filter(Ks > 1.140475) %>% 
  filter(Ks < 3.095325) %>% 
  select(X)

write.table(T.moresi_TEMO_WGD_paralogs, quote = F, file = "Tectaria_moresi_TEMO_WGD_paralogs_peak1.tsv")

#### Tectaria nayarii TENA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.nayarii_TENA <- read.delim("ks_distributions/Tectaria_nayarii.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.nayarii_TENA_filt <- T.nayarii_TENA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.nayarii_TENA_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.nayarii_TENA_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Tectaria nayarii (TENA)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.nayarii_TENA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.nayarii_TENA_filt.0 <- T.nayarii_TENA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.nayarii_TENA_normalmixEM <- normalmixEM(T.nayarii_TENA_filt.0$Ks, k=2)
# Fits to 1 WGD, comp 1 is overfit? 
summary(T.nayarii_TENA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.860324 0.139676
#mu     2.138817 0.166546
#sigma  0.976005 0.122665
#loglik at estimate:  -8948.52

# Generate list of paralogs 

T.nayarii_TENA_WGD_paralogs <- T.nayarii_TENA %>% 
  filter(Ks > 1.162812) %>% 
  filter(Ks < 3.114822) %>% 
  select(X)

write.table(T.nayarii_TENA_WGD_paralogs, quote = F, file = "Tectaria_nayarii_TENA_WGD_paralogs_peak1.tsv")

#### Tectaria subpedata TESU ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.subpedata_TESU <- read.delim("ks_distributions/Tectaria_subpedata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.subpedata_TESU_filt <- T.subpedata_TESU %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.subpedata_TESU_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.subpedata_TESU_filt, mapping = aes(x=Ks, ..scaled..*1100), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Tectaria subpedata (TESU)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.subpedata_TESU.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.subpedata_TESU_filt.0 <- T.subpedata_TESU_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.subpedata_TESU_normalmixEM <- normalmixEM(T.subpedata_TESU_filt.0$Ks, k=2)
# Fits to 1 WGD, comp 1 is overfit? 
summary(T.subpedata_TESU_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.0970654 0.902935
#mu     0.0876270 2.100936
#sigma  0.0795391 1.000635
#loglik at estimate:  -8160.315

# Generate list of paralogs 

T.subpedata_TESU_WGD_paralogs <- T.subpedata_TESU %>% 
  filter(Ks > 1.100301) %>% 
  filter(Ks < 3.101571) %>% 
  select(X)

write.table(T.subpedata_TESU_WGD_paralogs, quote = F, file = "Tectaria_subpedata_TESU_WGD_paralogs.tsv")

#### Thelypteris acuminata MROH ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.acuminata_MROH <- read.delim("ks_distributions/Thelypteris_acuminata_MROH.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.acuminata_MROH_filt <- T.acuminata_MROH %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.acuminata_MROH_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.acuminata_MROH_filt, mapping = aes(x=Ks, ..scaled..*600), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Thelypteris acuminata (MROH)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.acuminata_MROH.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.acuminata_MROH_filt.0 <- T.acuminata_MROH_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.acuminata_MROH_normalmixEM <- normalmixEM(T.acuminata_MROH_filt.0$Ks, k=2)
# Fits to 2 WGDs, ignore third component 
summary(T.acuminata_MROH_normalmixEM)
#summary of normalmixEM object:
#        comp 1   comp 2
#lambda 0.110481 0.889519
#mu     0.151772 2.150769
#sigma  0.113248 0.976703
#loglik at estimate:  -5867.427 

# Generate list of paralogs 

T.acuminata_MROH_WGD_paralogs <- T.acuminata_MROH %>% 
  filter(Ks > 1.174066) %>% 
  filter(Ks < 3.127472) %>% 
  select(X)

write.table(T.acuminata_MROH_WGD_paralogs, quote = F, file = "Thelypteris_acuminata_MROH_WGD_paralogs.tsv")

#### Thyrsopteris elegans EWXK ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.elegans_EWXK <- read.delim("ks_distributions/Thyrsopteris_elegans_EWXK.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.elegans_EWXK_filt <- T.elegans_EWXK %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.elegans_EWXK_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.elegans_EWXK_filt, mapping = aes(x=Ks, ..scaled..*1800), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Thyrsopteris elegans (EWXK)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.elegans_EWXK.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.elegans_EWXK_filt.0 <- T.elegans_EWXK_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.elegans_EWXK_normalmixEM <- normalmixEM(T.elegans_EWXK_filt.0$Ks, k=3)
# Fits to 1 WGD, ignore third component 
summary(T.elegans_EWXK_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.4604311 0.284134 0.255435
#mu     0.1440984 0.972687 2.617696
#sigma  0.0964686 0.461133 0.749865
#loglik at estimate:  -6274.888 

# Generate list of paralogs 

T.elegans_EWXK_WGD_paralogs <- T.elegans_EWXK %>% 
  filter(Ks > 0.511554) %>% 
  filter(Ks < 1.43382) %>% 
  select(X)

write.table(T.elegans_EWXK_WGD_paralogs, quote = F, file = "Thyrsopteris_elegans_EWXK_WGD_paralogs.tsv")

#### Tmesipteris parva ALVQ ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.parva_ALVQ <- read.delim("ks_distributions/Tmesipteris_parva_ALVQ.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.parva_ALVQ_filt <- T.parva_ALVQ %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.parva_ALVQ_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.parva_ALVQ_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Tmesipteris parva (ALVQ)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.parva_ALVQ.png", height = 5, width = 8, dpi = 300)

## Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.parva_ALVQ_filt.0 <- T.parva_ALVQ_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.parva_ALVQ_normalmixEM <- normalmixEM(T.parva_ALVQ_filt.0$Ks, k=3)
# Fits to 1 WGD (maybe 2??)
summary(T.parva_ALVQ_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.337998 0.277882 0.384120
#mu     0.161568 0.741274 2.366531
#sigma  0.116923 0.357751 0.851766
#loglik at estimate:  -7696.663

# Generate list of paralogs 

T.parva_ALVQ_WGD_paralogs <- T.parva_ALVQ %>% 
  filter(Ks > 0.383523) %>% 
  filter(Ks < 1.099025) %>% 
  select(X)

write.table(T.parva_ALVQ_WGD_paralogs, quote = F, file = "Tmesipteris_parva_ALVQ_WGD_paralogs.tsv")

#### Tmesipteris tannensis TMTA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.tannensis_TMTA <- read.delim("ks_distributions/Tmesipteris_tannensis.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.tannensis_TMTA_filt <- T.tannensis_TMTA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.tannensis_TMTA_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.tannensis_TMTA_filt, mapping = aes(x=Ks, ..scaled..*1200), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Tmesipteris tannensis (TMTA)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.tannensis_TMTA.png", height = 5, width = 8, dpi = 300)

## Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.tannensis_TMTA_filt.0 <- T.tannensis_TMTA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.tannensis_TMTA_normalmixEM <- normalmixEM(T.tannensis_TMTA_filt.0$Ks, k=3)
# Fits to 1 WGD (maybe 2??)
summary(T.tannensis_TMTA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2   comp 3
#lambda 0.360128 0.318386 0.321486
#mu     0.166684 0.835021 2.443975
#sigma  0.120761 0.405105 0.831017
#loglik at estimate:  -7837.455

# Generate list of paralogs 

T.tannensis_TMTA_WGD_paralogs <- T.tannensis_TMTA %>% 
  filter(Ks > 0.429916) %>% 
  filter(Ks < 1.240126) %>% 
  select(X)

write.table(T.tannensis_TMTA_WGD_paralogs, quote = F, file = "Tmesipteris_tannensis_TMTA_WGD_paralogs.tsv")

#### Trichomanes badium TRBA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

T.badium_TRBA <- read.delim("ks_distributions/Trichomanes_badium.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

T.badium_TRBA_filt <- T.badium_TRBA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(T.badium_TRBA_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(T.badium_TRBA_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Trichomanes badium (TRBA)") + theme(plot.title = element_text(face = "italic"))

ggsave("T.badium_TRBA.png", height = 5, width = 8, dpi = 300)

## Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

T.badium_TRBA_filt.0 <- T.badium_TRBA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

T.badium_TRBA_normalmixEM <- normalmixEM(T.badium_TRBA_filt.0$Ks, k=2, maxit = 2000)
# Fits to 1 WGD (maybe 2??)
summary(T.badium_TRBA_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.362337 0.637663
#mu     0.767475 2.458376
#sigma  0.482981 0.841804
#loglik at estimate:  -12787.11 

# Generate list of paralogs 

T.badium_TRBA_WGD_paralogs <- T.badium_TRBA %>% 
  filter(Ks > 1.616572) %>% 
  filter(Ks < 3.30018) %>% 
  select(X)

write.table(T.badium_TRBA_WGD_paralogs, quote = F, file = "Trichomanes_badium_TRBA_WGD_paralogs.tsv")

#### Vandenboschia striata VAST ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

V.striata_VAST <- read.delim("ks_distributions/Vandenboschia_striata.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

V.striata_VAST_filt <- V.striata_VAST %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(V.striata_VAST_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(V.striata_VAST_filt, mapping = aes(x=Ks, ..scaled..*500), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Vandenboschia striata (VAST)") + theme(plot.title = element_text(face = "italic"))

ggsave("V.striata_VAST.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD

#### Vittaria appalachiana NDUV ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

V.appalachiana_NDUV <- read.delim("ks_distributions/Vittaria_appalachiana_NDUV.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

V.appalachiana_NDUV_filt <- V.appalachiana_NDUV %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(V.appalachiana_NDUV_filt, mapping = aes(x=Ks), fill = "purple", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(V.appalachiana_NDUV_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "purple", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Vittaria appalachiana (NDUV)") + theme(plot.title = element_text(face = "italic"))

ggsave("V.appalachiana_NDUV.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

V.appalachiana_NDUV_filt.0 <- V.appalachiana_NDUV_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

V.appalachiana_NDUV_normalmixEM <- normalmixEM(V.appalachiana_NDUV_filt.0$Ks, k=2)
# Fits to 2 WGDs 
summary(V.appalachiana_NDUV_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.443437 0.556563
#mu     0.436252 2.422379
#sigma  0.287369 0.881995
#loglik at estimate:  -9297.802 

# Generate list of paralogs 

V.appalachiana_NDUV_WGD_paralogs_peak1 <- V.appalachiana_NDUV %>% 
  filter(Ks > 0.148883) %>% 
  filter(Ks < 0.723621) %>% 
  select(X)

write.table(V.appalachiana_NDUV_WGD_paralogs_peak1, quote = F, file = "Vittaria_appalachiana_NDUV_WGD_paralogs_peak1.tsv")

V.appalachiana_NDUV_WGD_paralogs_peak2 <- V.appalachiana_NDUV %>% 
  filter(Ks > 1.540384) %>% 
  filter(Ks < 3.304374) %>% 
  select(X)

write.table(V.appalachiana_NDUV_WGD_paralogs_peak2, quote = F, file = "Vittaria_appalachiana_NDUV_WGD_paralogs_peak2.tsv")

#### Vittaria lineata SKYV ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

V.lineata_SKYV <- read.delim("ks_distributions/Vittaria_lineata_SKYV.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

V.lineata_SKYV_filt <- V.lineata_SKYV %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(V.lineata_SKYV_filt, mapping = aes(x=Ks), fill = "pink", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(V.lineata_SKYV_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "pink", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Vittaria lineata (SKYV)") + theme(plot.title = element_text(face = "italic"))

ggsave("V.lineata_SKYV.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

V.lineata_SKYV_filt.0 <- V.lineata_SKYV_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

V.lineata_SKYV_normalmixEM <- normalmixEM(V.lineata_SKYV_filt.0$Ks, k=2)
# Fits to 2 WGDs 
summary(V.lineata_SKYV_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.466075 0.533925
#mu     0.482018 2.405027
#sigma  0.325006 0.882997
#loglik at estimate:  -7569.02 

# Generate list of paralogs 

V.lineata_SKYV_WGD_paralogs_peak1 <- V.lineata_SKYV %>% 
  filter(Ks > 0.157012) %>% 
  filter(Ks < 0.807024) %>% 
  select(X)

write.table(V.lineata_SKYV_WGD_paralogs_peak1, quote = F, file = "Vittaria_lineata_SKYV_WGD_paralogs_peak1.tsv")

V.lineata_SKYV_WGD_paralogs_peak2 <- V.lineata_SKYV %>% 
  filter(Ks > 1.52203) %>% 
  filter(Ks < 3.288024) %>% 
  select(X)

write.table(V.lineata_SKYV_WGD_paralogs_peak2, quote = F, file = "Vittaria_lineata_SKYV_WGD_paralogs_peak2.tsv")

#### Woodsia ilvensis WOIL ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.ilvensis_WOIL <- read.delim("ks_distributions/Woodsia_ilvensis_WOIL.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.ilvensis_WOIL_filt <- W.ilvensis_WOIL %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.ilvensis_WOIL_filt, mapping = aes(x=Ks), fill = "red", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.ilvensis_WOIL_filt, mapping = aes(x=Ks, ..scaled..*1500), fill = "red", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodsia ilvensis (WOIL)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.ilvensis_WOIL.png", height = 5, width = 8, dpi = 300)

# No evidence of WGD 

#### Woodsia ilvensis YQEC ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.ilvensis_YQEC <- read.delim("ks_distributions/Woodsia_ilvensis_YQEC.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.ilvensis_YQEC_filt <- W.ilvensis_YQEC %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.ilvensis_YQEC_filt, mapping = aes(x=Ks), fill = "orange", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.ilvensis_YQEC_filt, mapping = aes(x=Ks, ..scaled..*900), fill = "orange", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodsia ilvensis (YQEC)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.ilvensis_YQEC.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

W.ilvensis_YQEC_filt.0 <- W.ilvensis_YQEC_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

W.ilvensis_YQEC_normalmixEM <- normalmixEM(W.ilvensis_YQEC_filt.0$Ks, k=2)
# Fits to 1 WGD, ignore second comp
summary(W.ilvensis_YQEC_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.833654 0.166346
#mu     1.561493 3.403225
#sigma  0.947775 0.365777
#loglik at estimate:  -6007.91 

# Generate list of paralogs 

W.ilvensis_YQEC_WGD_paralogs <- W.ilvensis_YQEC %>% 
  filter(Ks > 0.613718) %>% 
  filter(Ks < 2.509268) %>% 
  select(X)

write.table(W.ilvensis_YQEC_WGD_paralogs, quote = F, file = "Woodsia_ilvensis_YQEC_WGD_paralogs.tsv")

#### Woodsia polystichoides WOPO ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.polysitchoides_WOPO <- read.delim("ks_distributions/Woodsia_polystichoides.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.polysitchoides_WOPO_filt <- W.polysitchoides_WOPO %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.polysitchoides_WOPO_filt, mapping = aes(x=Ks), fill = "yellow", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.polysitchoides_WOPO_filt, mapping = aes(x=Ks, ..scaled..*300), fill = "yellow", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodsia polysitchoides (WOPO)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.polysitchoides_WOPO.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

W.polysitchoides_WOPO_filt.0 <- W.polysitchoides_WOPO_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

W.polysitchoides_WOPO_normalmixEM <- normalmixEM(W.polysitchoides_WOPO_filt.0$Ks, k=2)
# Fits to 1 WGD 
summary(W.polysitchoides_WOPO_normalmixEM)
#summary of normalmixEM object:
#  comp 1   comp 2
#lambda 0.0961639 0.903836
#mu     0.2346062 2.185603
#sigma  0.1578547 0.934814
#loglik at estimate:  -6391.354 

# Generate list of paralogs 

W.polysitchoides_WOPO_WGD_paralogs <- W.polysitchoides_WOPO %>% 
  filter(Ks > 1.250789) %>% 
  filter(Ks < 3.120417) %>% 
  select(X)

write.table(W.polysitchoides_WOPO_WGD_paralogs, quote = F, file = "Woodsia_polysitchoides_WOPO_WGD_paralogs.tsv")

#### Woodwardia japonica WOJA ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.japonica_WOJA <- read.delim("ks_distributions/Woodwardia_japonica.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.japonica_WOJA_filt <- W.japonica_WOJA %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.japonica_WOJA_filt, mapping = aes(x=Ks), fill = "cyan", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.japonica_WOJA_filt, mapping = aes(x=Ks, ..scaled..*700), fill = "cyan", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodwardia japonica (WOJA)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.japonica_WOJA.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

W.japonica_WOJA_filt.0 <- W.japonica_WOJA_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

W.japonica_WOJA_normalmixEM <- normalmixEM(W.japonica_WOJA_filt.0$Ks, k=2)
# Fits to 1 WGD 
summary(W.japonica_WOJA_normalmixEM)
#summary of normalmixEM object:
#       comp 1  comp 2
#lambda 0.165790 0.83421
#mu     0.157092 2.08278
#sigma  0.123759 1.00574
#loglik at estimate:  -9571.334  

# Generate list of paralogs 

W.japonica_WOJA_WGD_paralogs <- W.japonica_WOJA %>% 
  filter(Ks > 1.07704) %>% 
  filter(Ks < 3.08852) %>% 
  select(X)

write.table(W.japonica_WOJA_WGD_paralogs, quote = F, file = "Woodwardia_japonica_WOJA_WGD_paralogs.tsv")

#### Woodwaria prolifera WOPR ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.prolifera_WOPR <- read.delim("ks_distributions/Woodwardia_prolifera.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.prolifera_WOPR_filt <- W.prolifera_WOPR %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.prolifera_WOPR_filt, mapping = aes(x=Ks), fill = "blue", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.prolifera_WOPR_filt, mapping = aes(x=Ks, ..scaled..*1000), fill = "blue", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodwardia prolifera (WOPR)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.prolifera_WOPR.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

W.prolifera_WOPR_filt.0 <- W.prolifera_WOPR_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

W.prolifera_WOPR_normalmixEM <- normalmixEM(W.prolifera_WOPR_filt.0$Ks, k=2)
# Fits to 1 WGD 
summary(W.prolifera_WOPR_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.145996 0.854004
#mu     0.124887 2.100409
#sigma  0.117699 0.979222
#loglik at estimate:  -10506.92  

# Generate list of paralogs 

W.prolifera_WOPR_WGD_paralogs <- W.prolifera_WOPR %>% 
  filter(Ks > 1.21187) %>% 
  filter(Ks < 3.079631) %>% 
  select(X)

write.table(W.prolifera_WOPR_WGD_paralogs, quote = F, file = "Woodwardia_prolifera_WOPR_WGD_paralogs.tsv")

#### Woodsia scopulina YJJY ####

# Family-wise distributions from wgd (Zwaenepoel and Van de Peer 2019; https://github.com/arzwa/wgd)

W.scopulina_YJJY <- read.delim("ks_distributions/Woodsia_scopulina_YJJY.cds.ks.tsv") 

# Filter out saturated duplicates (Ks < 4)

W.scopulina_YJJY_filt <- W.scopulina_YJJY %>% 
  select(Ks) %>%  
  filter(Ks < 4)

ggplot() + geom_histogram(W.scopulina_YJJY_filt, mapping = aes(x=Ks), fill = "darkgreen", color ="black", alpha = 0.7, binwidth = 0.1) + 
  geom_density(W.scopulina_YJJY_filt, mapping = aes(x=Ks, ..scaled..*800), fill = "darkgreen", alpha = 0.4, adjust=0.6) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.001, name = "Density")) +
  xlim(-0.1, 4) + theme_classic() + ylab("Count") + ggtitle("Woodsia scopulina (YJJY)") + theme(plot.title = element_text(face = "italic"))

ggsave("W.scopulina_YJJY.png", height = 5, width = 8, dpi = 300)

# Run mixture model with mixtools (Benaglia et al. 2009; https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
# Re-filter for recent tandem duplicates 

W.scopulina_YJJY_filt.0 <- W.scopulina_YJJY_filt %>% 
  select(Ks) %>%
  filter(Ks > 0)

W.scopulina_YJJY_normalmixEM <- normalmixEM(W.scopulina_YJJY_filt.0$Ks, k=2)
# Fits to 1 WGD, ignore second comp
summary(W.scopulina_YJJY_normalmixEM)
#summary of normalmixEM object:
#         comp 1   comp 2
#lambda 0.120299 0.879701
#mu     0.125141 2.180207
#sigma  0.105676 0.960208
#loglik at estimate:  -6008.843 

# Generate list of paralogs 

W.scopulina_YJJY_WGD_paralogs <- W.scopulina_YJJY %>% 
  filter(Ks > 1.219999) %>% 
  filter(Ks < 3.140415) %>% 
  select(X)

write.table(W.scopulina_YJJY_WGD_paralogs, quote = F, file = "Woodsia_scopulina_YJJY_WGD_paralogs.tsv")
