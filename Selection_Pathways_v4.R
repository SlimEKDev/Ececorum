#!/usr/bin/env Rscript

setwd("~/Ececorum/datastxt4")

# Signale la possibilité d'avoir des arguments et les transforme en numériques
args = as.numeric(commandArgs(TRUE))

# Lecture des heatmaps téléchargées sous forme de tableau tabulés et concaténation sous un seul fichier
# Les warnings sont dûs :
# 1- à la création de out.file vide => ajout de NA auto 
# 2- et à la lecture des EC avec une virgule dans le nom qui entrainent un décalage, celui-ci étant corrigé par la lecture du .txt

# Création d'un fichier vide pour éviter l'erreur à la première lecture
out.file <- ""

# Création d'une liste contenant tous les noms des fichiers .txt
file_names <- dir(path = ".", pattern =".txt")

# Lecture des .txt, réarrangement des données et concaténation
for(i in 1:length(file_names)) {
  file2 <- read.csv(file = file_names[i], header=TRUE, row.names=NULL, sep="\t" , fill = TRUE) # Lecture du fichier
  file2$Pathway.ID <- file_names[i] # ajout des noms des pathways
  file2 <- file2[,c(length(file2),1:length(file2)-1)] # déplacement du nom du pathway en 2ème colonne
  out.file <- rbind(out.file, file2) #concaténation
}

# Elimine la 1ère ligne créée par le fichier vide
out.file <- out.file[2:dim(out.file)[1],]

# Elimine les colonnes correspondant aux souches d'E. ceco redondante pour la souche DSM20682
out.file <- out.file[,-c(16:18)]

# Renommer les Pathways correctement, sans le .txt
vec = c()
for(i in 1:length(out.file$Pathway.ID)){
  vec <- c(vec, substr(out.file$Pathway.ID[i],1,nchar(out.file$Pathway.ID[i])-4))
}
out.file$Pathway.ID <- vec

# Détermination des colonnes de cecorum
ColCeco <- grep('cecorum', colnames(out.file))

# Détermination des colonnes de columbae
ColColum <- grep('columbae', colnames(out.file))

# Détermination de la taille totale de l'échantillon de génomes
ColEnter <- grep('Enterococcus', colnames(out.file))

# Construction d'une df uniquement pour E. cecorum
CecoTab <- out.file[,c(1, 2, ColCeco)]

# Comptage des génomes possédant au moins 1 exemplaire de l'E.C.
CecoTab$Nb.Hit <- apply(CecoTab[,3:length(CecoTab)],1, function(x) sum(x>=1))

# Attente d'un argument sinon le réclame
if (is.na(args[1]) == TRUE) {
  TolCeco <- as.numeric(readline(prompt=paste("Combien de E. cecorum doivent posséder le pathway sur", length(ColCeco),"? ", sep=" "))) #modulateur de tolérance des E. ceco
} else {
  TolCeco <- args[1]
}

# Sélection des EC présentss chez E. cecorum (ou une fraction selon l'argument)
CecoTab <- CecoTab[CecoTab$Nb.Hit >= TolCeco, ]

# Construction d'une df uniquement pour E. columbae
ColTab <- out.file[,c(1, 2, ColColum)]

# Comptage des génomes possédant au moins 1 exemplaire de l'E.C.
ColTab$Nb.Hit <- apply(ColTab[,3:length(ColTab)],1, function(x) sum(x>=1))

# Attente d'un argument (modulation de la stringence) sinon le réclame
if (is.na(args[2]) == TRUE) {
  TamisCol <- as.numeric(readline(prompt=paste("Combien de E. columbae peuvent posséder le pathway sur", length(ColColum),"? ", sep=" "))) #modulateur de tolérance des E. colum
} else {
  TamisCol <- args[2]
}

# Sélection des EC absents chez columbae (ou une fraction selon l'argument)
ColTab <- ColTab[ColTab$Nb.Hit <= TamisCol, ]

# Construction d'une df avec tous les Enterococcus sans cecorum et columbae
SppTab <- out.file[,-c(ColCeco, ColColum)]

# Comptage des génomes possédant au moins 1 exemplaire de l'E.C.
SppTab$Nb.Hit <- apply(SppTab[,3:length(SppTab)],1, function(x) sum(x>=1))

# Détermination des Enterococcus spp sans ceco et sans columbae
ColSpp <- (length(ColEnter)-length(ColCeco)-length(ColColum))

# Attente d'un argument sinon le réclame
if (is.na(args[3]) == TRUE) {
  TamisSpp <- as.numeric(readline(prompt=paste("Combien de E. spp peuvent posséder le pathway sur",ColSpp,"? ", sep = " "))) #modulateur de tolérance des E. spp
} else {
  TamisSpp <- args[3]
}

# Sélection des EC absents chez les E. spp (ou une fraction selon l'argument)
SppTab <- SppTab[SppTab$Nb.Hit <= TamisSpp, ]

# Application des 3 tamis
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library") #librairie pour utiliser semi_join

# Tamis pour E. cecorum
out.file <- semi_join(out.file,CecoTab, by="Protein.Families.Genomes")
CecoTabtemp <- semi_join(CecoTab, out.file, by="Protein.Families.Genomes")
out.file$Nb.ECecorum <- CecoTabtemp$Nb.Hit

# Tamis pour E. columbae
out.file <- semi_join(out.file,ColTab, by="Protein.Families.Genomes")
ColTabtemp <- semi_join(ColTab, out.file, by="Protein.Families.Genomes")
out.file$Nb.EColumbae <- ColTabtemp$Nb.Hit

# Tamis pour E. spp
out.file <- semi_join(out.file,SppTab, by="Protein.Families.Genomes")
SppTabtemp <- semi_join(SppTab, out.file, by="Protein.Families.Genomes")
out.file$Nb.ESpp <- SppTabtemp$Nb.Hit

# Ajout du nom du pathway
PathEcec <- read.csv(file = "../Ececorum_all17_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE, colClasses=c(rep("factor",5))) # Lecture du .csv contenant les noms des Pathways de E. cecorum
CatPathID <- PathEcec[,c(1,3)] # Sélection des colonnes contenant l'ID du pathway et celle contenant son nom
CatPathID$Pathway.ID <- sapply(CatPathID$Pathway.ID, function(x) toString(x)) # Transformation des ID de pathways pour les rendre compatibles avec ceux d'out.file
out.file <- merge(out.file,CatPathID, by="Pathway.ID") # Ajout de la df avec le nom du pathway à notre out.file

# Création d'une colonne d'interrogation de KEGG
out.file$KEGG.link <- sapply(out.file$Protein.Families.Genomes, function(x) paste("http://www.genome.jp/dbget-bin/www_bget?ec:",gsub(" .*", "", x),sep=""))

# Réarrangement des colonnes
out.file <- out.file[,c(1,length(out.file)-1,2,length(out.file),length(out.file)-4,length(out.file)-3,length(out.file)-2,3:(length(out.file)-5))]

#Naming avec parametres
name_outfile <- paste("PathwaysOfInterest-",nrow(out.file),"_","Ececo-",round((TolCeco/length(ColCeco))*100,2),"%_ECol-",round((TamisCol/length(ColColum))*100,2),"%_TamisSpp-",round((TamisSpp/(length(ColEnter)-length(ColCeco)-length(ColColum)))*100,2),"%.csv", sep = "")
path <- paste("./Resultats4/",name_outfile, sep = "")

# Ecriture du .csv de sortie
write.csv(out.file, file = path)

# Signal de la génération du fichier
cat("Fichier enregistré :",name_outfile)
