setwd("~/Ececorum")

# Chargement souches complètes cecorum
PathEcecCplt <- read.csv(file = "Ececorum_gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)

CandidPathEcecCplt <- PathEcecCplt[PathEcecCplt$Unique.Genome.Count == 6 , ]
# Chargement Comparaison Pathway pour E. spp
PathEspp41 <- read.csv(file = "Espp_41-gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
# Chargement Comparaison Pathway pour E. spp sans E. cecorum
PathEspp41mEceco <- read.csv(file = "Espp_41-gnm-cplt_mEceco_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
# Séléction Candidat si Nb EC augmente
CandidPathEspp41 <- PathEspp41[ PathEspp41$Unique.EC.Count > PathEspp41mEceco$Unique.EC.Count , ]
# Sélection si présence dans pathway candidat de E cecorum
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
PathECSup <- semi_join(CandidPathEspp41, CandidPathEcecCplt, by = "Pathway.ID")

CandidPathEcecCplt <- semi_join(CandidPathEcecCplt, PathECSup, by = "Pathway.ID")
PathECSup$Prev.EC.Count <- CandidPathEcecCplt$Unique.EC.Count
CCC <- PathECSup[,c(1,2,3,10,7)]
#BBB <- data.frame(v=1:4, ch=c(PathECSup$Pathway.ID, PathECSup$Pathway.Name, PathECSup$Pathway.Class, PathECSup$Prev.EC.Count, PathECSup$Unique.EC.Count))
#AAA <- anti_join(CandidPathEspp41, PathECSup)
write.csv(PathECSup, file = "Pathw-Candid_Cribl-Sup-EC.csv")
