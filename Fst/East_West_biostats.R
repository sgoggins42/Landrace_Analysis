source("hierfstat.R")

#Import East and West biostatistics
EW_basicstats <- read.csv("~/Desktop/EW_basicstats.csv")
EW_Fst <- EW_basicstats[,c("X","Fst")]
#Select SNP location info
SNP_by_chro <- GeneticMap_iSelect_9k[,c("SNP","cm","Chromosome")]
#merge location of SNP to FST values
SNP_FST.df <- merge(x = EW_Fst, y = SNP_by_chro, by.x = "X", by.y = "SNP")

#import inversed and edited data
tgeno <- read.csv("~/Desktop/Landrace_Analysis/Worked_Datasets/tgeno.csv")
#merge
outliers.in.samples <- merge(x = Outliers_FST, y = tgeno, by.x = "X", by.y = "SNP.ID")

#Identify Fst Outliers by 95%
Fst_values <- SNP_FST.df$Fst
Fst_range <- range(Fst_values,na.rm=TRUE) * .95
Fst_range <- as.data.frame(Fst_range)
Outlier_Fst <- subset(SNP_FST.df, Fst > Fst_range[2,1])
outliers.in.samples <- merge(x = Outlier_Fst, y = tgeno, by.x = "X", by.y = "SNP.ID")

#Visualize FST Distribution by Chromosome
with(SNP_FST.df, plot(x = jitter(Chromosome,2),y = Fst,
                      cex = .75, col = ifelse(Fst > Fst_range[2,1], "red", "black"),
                      xlab = "Chromosome", ylab = "Fst value", main = "Outlier Fst Values"))
abline(h = Fst_range[2,1],lty = 6, lwd = .5, col = "red")

#map of chromosome 2 outliers
library("maps")
library("maptools")

quartz()
data("wrld_simpl")
plot(wrld_simpl,xlim=c(-30,140),ylim=c(-8,64),col="lemonchiffon2",main="Barley Landraces")
box()

locgenoEW.df$side <- ifelse(test = locgenoEW.df$X11_10317 == 1, yes = "1", no = "2")
SNP.ID.2 <- locgenoEW.df[,c("Accession.ID", "side")]
loc.df <- loc.df[,c(1:5)]
SNP.ID.2 <- merge(x = SNP.ID.2, y = loc.df, by.x = "Accession.ID", by.y = "Accession.ID")
with(SNP.ID.2, points(x = Longitude,y = Latitude, cex = .75,
                      col = ifelse(side == "1", "blue","red")))
abline(v=48) #Zagaros 
minor_allele <- table(locgenoEW.df$side)

legend("bottomleft", bg = "transparent",
       c("Zagaros Mountains","Minor allele", "Major allele"),
       col=c("black","blue","red"), lwd=1, lty=c(1,NA,NA),
       pch = c(NA,21,21),cex = 0.5)
