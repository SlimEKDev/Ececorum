setwd("~/Ececorum/datastxt")

#Lecture de heatmap DL sous forme de tableau tabulés
# Les warnings sont dûs à la création de out.file vide => ajout de NA auto 
# et à la lecture des EC avec une virgule dans le nom qui entrainent un décalage, corrigé par la lecture du .txt
out.file <- ""
file_names <- dir(path = ".", pattern =".txt")
for(i in 1:length(file_names)) {
  file2 <- read.csv(file = file_names[i], header=TRUE, row.names=NULL, sep="\t" , fill = TRUE)
  file2$Pathway.ID <- file_names[i] # ajout des noms des pathways
  file2 <- file2[,c(58,1:57)] # déplacement du nom du pathway en 1ère colonne
  out.file <- rbind(out.file, file2)
}

# Eliminer la 1ère ligne créée par le fichier vide
out.file=out.file[2:dim(out.file)[1],]

# Renommer les Pathways correctement, sans le .txt
vec = c()
for(i in 1:length(out.file$Pathway.ID)){
  vec <- c(vec, substr(out.file$Pathway.ID[i],1,nchar(out.file$Pathway.ID[i])-4))
}
out.file$Pathway.ID <- vec

# Détermination des colonnes de cecorum et de columbae
ColCeco <- grep('cecorum', colnames(out.file))
ColColum <- grep('columbae', colnames(out.file))

# Construction d'une df uniquement pour E. ceco
Ceco1 <- out.file[,c(1, 2, ColCeco)]
# Sélection des EC que possèdent tous les E. ceco (ou une grande fraction)
Ceco1$Nb.Zero <- apply(Ceco1,1, function(x) sum(x==0))
TolCeco <- 0/17 #modulateur de tolérance des E. ceco
Ceco1 <- Ceco1[Ceco1$Nb.Zero <= (length(Ceco1)-3)*TolCeco, ]

# Construction d'une df uniquement pour E. columbae
Ceco2 <- out.file[,c(1, 2, ColColum)]
# Sélection des EC absents chez columbae
Ceco2$Nb.Zero <- apply(Ceco2,1, function(x) sum(x==0))
TamisCol <- 4/4 #modulateur de stringence du nombre de E. columbae
Ceco2 <- Ceco2[Ceco2$Nb.Zero >= (length(Ceco2)-3)*TamisCol, ]

# Construction d'une df avec tous les Enterococcus sans cecorum et columbae
Ceco3 <- out.file[,-c(ColCeco, ColColum)]
# Sélection des EC absents chez tous les E. spp
Ceco3$Nb.Zero <- apply(Ceco3,1, function(x) sum(x==0))
TamisSpp <- 1 #modulateur de stringence du nombre E. spp (cecorum et columbae exclus)
Ceco3 <- Ceco3[Ceco3$Nb.Zero >= (length(Ceco3)-3)*TamisSpp, ]

# Application des 3 tamis
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
out.file <- semi_join(out.file,Ceco1)
out.file <- semi_join(out.file,Ceco2)
out.file <- semi_join(out.file,Ceco3)

# Ajout des noms des Pathways
PathEcec <- read.csv(file = "../Ececorum_all17_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
CatPathID <- PathEcec[,c(1,3)]
temp <- sapply(CatPathID$Pathway.ID, function(x) toString(x))
CatPathID$Pathway.ID <- temp
out.file <- inner_join(out.file,CatPathID)
out.file <- out.file[,c(1,59,2:58)] # déplacement du nom du pathway en 2ème colonne

# Naming avec paramètres
name_outfile <- paste("PathwaysOfInterest_Ececo-",round((1-TolCeco)*100,2),"%_ECol-",round((1-TamisCol)*100,2),"%_TamisSpp-",round((1-TamisSpp)*100,2),"%.csv", sep = "")
write.csv(out.file, file = name_outfile)
