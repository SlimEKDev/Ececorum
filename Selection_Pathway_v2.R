setwd("~/Ececorum/datastxt")

# Lecture de heatmap DL sous forme de tableau tabulés
# Les warnings sont dûs à la création de out.file vide => ajout de NA auto 
# et à la lecture des EC avec une virgule dans le nom qui entrainent un décalage, corrigé par la lecture du .txt
out.file <- ""
file_names <- dir(path = ".", pattern =".txt")
for(i in 1:length(file_names)) {
  file2 <- read.csv(file = file_names[i], header=TRUE, row.names=NULL, sep="\t" , fill = TRUE)
  file2$Pathway.ID <- file_names[i] # ajout des noms des pathways
  file2 <- file2[,c(length(file2),1:length(file2)-1)] # déplacement du nom du pathway en 2ème colonne
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

#détermination des colonnes de cecorum et de columbae
ColCeco <- grep('cecorum', colnames(out.file))
ColColum <- grep('columbae', colnames(out.file))

#détermination de la taille de l'échantillon de souches
ColEnter <- grep('Enterococcus', colnames(out.file))

#Construction d'une df uniquement pour E. ceco
Ceco1 <- out.file[,c(1, 2, ColCeco)]
#sélection des EC que possèdent tous les E. ceco (ou une grande fraction)
Ceco1$Nb.Hit <- apply(Ceco1[,3:length(Ceco1)],1, function(x) sum(x>=1))
TolCeco <- as.numeric(readline(prompt=paste("Combien de E. cecorum doivent posséder le pathway sur", length(ColCeco),"? ", sep=" "))) #modulateur de tolérance des E. ceco
Ceco1 <- Ceco1[Ceco1$Nb.Hit >= TolCeco, ]
#(length(Ceco1)-3) nombre de cecorum

#Construction d'une df uniquement pour E. columbae
Ceco2 <- out.file[,c(1, 2, ColColum)]
#sélection des EC absents chez columbae
Ceco2$Nb.Hit <- apply(Ceco2[,3:length(Ceco2)],1, function(x) sum(x>=1))
TamisCol <- as.numeric(readline(prompt=paste("Combien de E. columbae peuvent posséder le pathway sur", length(ColColum),"? ", sep=" "))) #modulateur de tolérance des E. colum
#TamisCol <- 0/4 #modulateur de stringence du nombre de E. columbae
Ceco2 <- Ceco2[Ceco2$Nb.Hit <= TamisCol, ]

#Construction d'une df avec tous les Enterococcus sans cecorum et columbae
Ceco3 <- out.file[,-c(ColCeco, ColColum)]
#sélection des EC absents chez tous les E. spp
Ceco3$Nb.Hit <- apply(Ceco3[,3:length(Ceco3)],1, function(x) sum(x>=1))
#détermination des Enterococcus spp sans ceco et sans columbae
ColSpp <- (length(ColEnter)-length(ColCeco)-length(ColColum))
TamisSpp <- as.numeric(readline(prompt=paste("Combien de E. spp peuvent posséder le pathway sur",ColSpp,"? ", sep = " "))) #modulateur de tolérance des E. spp
#TamisSpp <- 0 #modulateur de stringence du nombre E. spp (cecorum et columbae exclus)
Ceco3 <- Ceco3[Ceco3$Nb.Hit <= TamisSpp, ]

#application des 3 tamis
library("dplyr", lib.loc="/usr/local/public/R-3.2.0/lib64/R/library")
out.file <- semi_join(out.file,Ceco1)
out.file$Nb.ECecorum <- Ceco1$Nb.Hit
out.file <- semi_join(out.file,Ceco2)
Ceco2temp <- semi_join(Ceco2, out.file)
out.file$Nb.EColumbae <- Ceco2temp$Nb.Hit
out.file <- semi_join(out.file,Ceco3)
Ceco3temp <- semi_join(Ceco3, out.file)
out.file$Nb.ESpp <- Ceco3temp$Nb.Hit

PathEcec <- read.csv(file = "../Ececorum_all17_PATRIC_pathways.csv", header = TRUE, sep = "," , fill = TRUE)
CatPathID <- PathEcec[,c(1,3)]
temp <- sapply(CatPathID$Pathway.ID, function(x) toString(x)) 
CatPathID$Pathway.ID <- temp
out.file <- inner_join(out.file,CatPathID)
#out.file <- out.file[,c(1,62,59,60,61,3:58)] # déplacement du nom du pathway en 1ère colonne
out.file <- out.file[,c(1,length(out.file),length(out.file)-3,length(out.file)-2,length(out.file)-1,3:(length(out.file)-4))] # déplacement du nom du pathway en 1ère colonne

#Naming avec parametres
name_outfile <- paste("PathwaysOfInterest_Ececo-",round((TolCeco/length(ColCeco))*100,2),"%_ECol-",round((TamisCol/length(ColColum))*100,2),"%_TamisSpp-",round((TamisSpp/(length(ColEnter)-length(ColCeco)-length(ColColum)))*100,2),"%.csv", sep = "")
write.csv(out.file, file = name_outfile)
cat("Fichier enregistré :",name_outfile) 
