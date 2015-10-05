# source("hierfstat.R")

#Import East and West biostatistics
EW_basicstats <- read.csv("~/Desktop/EW_basicstats.csv")
EW_Fst <- EW_basicstats[,c("X","Fst")]
#Select SNP location info
SNP_by_chro <- GeneticMap_iSelect_9k[,c("SNP","cm","Chromosome")]
#merge location of SNP to FST values
SNP_FST.df <- merge(x = EW_Fst, y = SNP_by_chro, by.x = "X", by.y = "SNP")
#Visualize Fst values
with(SNP_FST.df, plot(x = jitter(Chromosome,2),y = Fst,
                      cex = .75, col = ifelse(Fst > .5, "red", "black"), pch = ifelse(Fst > .5, 19, 21),
                      xlab = "Chromosome", ylab = "Fst value", main = "Outlier Fst Values"))
abline(h=.5, lty = 6, lwd = .5, col = "red")
#table of outliers
Outliers <- subset(SNP_FST.df, Fst > .5)
#import inversed and edited data
tgeno <- read.csv("~/Desktop/Landrace_Analysis/Worked_Datasets/tgeno.csv")
#merge
outliers.in.samples <- merge(x = Outliers, y = tgeno, by.x = "X", by.y = "SNP.ID")
#2H pericentric inversion: ~67-74cm
chro.2.inversion <- subset(Outliers, Chromosome == 2)


#map of chromosome 2 outliers
library("maps")
library("maptools")

quartz()
data("wrld_simpl")
plot(wrld_simpl,xlim=c(-30,140),ylim=c(-8,64),col="lemonchiffon2",main="Barley Landraces")
box()

locgenoEW.df$side <- ifelse(test = locgenoEW.df$SCRI_RS_185513 == 2, yes = "1", no = "2")
SNP.ID.2 <- locgenoEW.df[,c("Accession.ID", "side")]
loc.df <- loc.df[,c(1:5)]
SNP.ID.2 <- merge(x = SNP.ID.2, y = loc.df, by.x = "Accession.ID", by.y = "Accession.ID")
with(SNP.ID.2, points(x = Longitude,y = Latitude, cex = .75,
                      col = ifelse(side == "1", "red","blue")))
abline(v=48) #Zagaros 

legend("bottomleft", bg = "transparent",
       c("Zagaros Mountains","Derived Chromosome 2", "Ancestral Chromosome 2"),
       col=c("black","blue","red"), lwd=1, lty=c(1,NA,NA),
       pch = c(NA,21,21),cex = 0.8)
