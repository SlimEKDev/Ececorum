setwd("~/Ececorum")
# Chargement souches complètes cecorum
PathEcecCplt <- read.csv(file = "Ececorum_gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
CandidPathEcecCplt <- PathEcecCplt[PathEcecCplt$Unique.Genome.Count == 6 & PathEcecCplt$EC.Conservation == 100, ]
# Chargement Comparaison Pathway pour E. spp
PathEspp41 <- read.csv(file = "Espp_41-gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
# Chargement Comparaison Pathway pour E. spp = pan-pathwome de E. spp, E. cecorum exclu
PathEspp41mEceco <- read.csv(file = "Espp_41-gnm-cplt_mEceco_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
# Séléction Candidat si EC.Conservation augmente
CandidPathEspp41 <- PathEspp41[PathEspp41$EC.Conservation >= PathEspp41mEceco$EC.Conservation , ]
summary(CandidPathEcecCplt)
summary(CandidPathEspp41)
# Sélection si présence dans pathway candidat de E cecorum
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
PathECAmelio <- semi_join(CandidPathEspp41, CandidPathEcecCplt, by = 'Pathway.ID')
write.csv(PathECAmelio, file = "Pathw-Candid_Cribl-Amelio-EC.csv")
