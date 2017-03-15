setwd("~/Ececorum")
# Chargement Comparaison Pathway E. cecorum
GenEcec <- read.csv(file = "Ececorum_all_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEcec)
# Selection Pathway si retrouvé 16 fois (soit tous les génomes sélectionnés) = core-métabolome de E. cecorum
CoreGenEcec <- GenEcec[ GenEcec[,5] == 16, c(1,2,5) ]
summary(CoreGenEcec)
# Chargement Comparaison Pathway pour E. spp
GenEspp <- read.csv(file = "Espp_gnm-cplt_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEspp)
# Sélection Pathway si retrouvé 37 fois (soit tous les génomes sélectionnés) = core-metabolome de E. spp
CoreGenEspp <- GenEspp[ GenEspp[,5] == 37, c(1,2,5) ]
summary(CoreGenEspp)
# Sélection des pathways uniques à E. cecorum en éliminant les redoncances entre les 2 tables
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
GenSpeEcec <- anti_join(CoreGenEcec,CoreGenEspp, by = "Pathway.ID")
# Sauvegarde de la table en fichier .csv
write.csv(GenSpeEcec, file = "Genome-specifique_E-cecorum.csv")

setwd("~/Ececorum")
# Chargement Comparaison Pathway E. cecorum
GenEcec <- read.csv(file = "Ececorum_all_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEcec)
#selection Pathway si retrouvé 16 fois (soit tous les génomes sélectionnés)= = core-métabolome de E. cecorum
CoreGenEcec <- GenEcec[ GenEcec[,5] == 16, c(1,2,5) ]
summary(CoreGenEcec)
# Chargement Comparaison Pathway pour E. spp = pan-metabolome E. spp, E. cecorum exclu
GenEsppmEceco <- read.csv(file = "Espp_gnm-cplt_mEceco_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
summary(GenEsppmEceco)
# Sélection colonnes 1,2,5 dans la table pan-metabolomede E. spp, E. cecorum exclu
PanGenEsppmEceco <- GenEsppmEceco[ , c(1,2,5) ]
summary(CoreGenEspp)
# Prise en compte uniquement des gènes variables 
VarGenEsppmEceco <- GenEsppmEceco[ GenEsppmEceco[,5] < 31, c(1,2,5) ]
# Détermination du nombre d'occurence des pathways "spécifiques" de E. cecorum dans le génome variable de E. spp
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
CommonGen <- semi_join(VarGenEsppmEceco, CoreGenEcec, by="Pathway.ID")
# Ecriture résultat dans un fichier .csv
write.csv(CommonGen, file = "NB-occ_Genome-unique_E-cecorum.csv")
