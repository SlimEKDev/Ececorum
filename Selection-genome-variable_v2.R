setwd("~/Ececorum")
# Chargement Comparaison Pathway E. cecorum
GenEcec <- read.csv(file = "Ececorum_all_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEcec)
# Sélection Pathway si retrouvé 16 fois (soit tous les génomes sélectionnés) = core-pathwome de E. cecorum
CoreGenEcec <- GenEcec[ GenEcec[,5] == 16, c(1,2,5,8) ]
summary(CoreGenEcec)

# Chargement Comparaison Pathway pour E. spp
GenEspp <- read.csv(file = "Espp_41-gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEspp)
# Sélection Pathway si retrouvé 41 fois (soit tous les génomes sélectionnés) = core-pathwome de E. spp
CoreGenEspp <- GenEspp[ GenEspp[,5] == 41, c(1,2,5,8) ]
summary(CoreGenEspp)
# Sélection des pathways uniques à E. cecorum en éliminant les redondances entre les 2 tables
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
GenSpeEcec <- anti_join(CoreGenEcec,CoreGenEspp, by = "Pathway.ID")
# Sauvegarde de la table en fichier .csv
write.csv(GenSpeEcec, file = "Genome-specifique_E-cecorum2.csv")

# Chargement Comparaison Pathway pour E. spp = pan-pathwome de E. spp, E. cecorum exclu
GenEsppmEceco <- read.csv(file = "Espp_41-gnm-cplt_mEceco_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEsppmEceco)
# Sélection colonnes 1,2,5,8 dans la table pan-pathwome de E. spp, E. cecorum exclu
PanGenEsppmEceco <- GenEsppmEceco[ , c(1,2,5,8) ]
summary(PanGenEsppmEceco)
# Prise en compte uniquement des gènes variables 
VarGenEsppmEceco <- GenEsppmEceco[ GenEsppmEceco[,5] < 35, c(1,2,5,8) ]
summary(VarGenEsppmEceco)
CommonGen <- semi_join(VarGenEsppmEceco, CoreGenEcec, by="Pathway.ID")
# Ecriture résultat dans un fichier .csv
write.csv(CommonGen, file = "NB-occ_Pathwome-unique_E-cecorum2.csv")
