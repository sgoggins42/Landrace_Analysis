#install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import dataset
`803_landraces_KML` <- read.csv("~/Desktop/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")
Land_6152_SNPs_AB <- read.delim("~/Desktop/Landrace_Analysis/Worked_Datasets/Land_6152_SNPs_AB.txt")
GeneticMap_iSelect_9k <- read.delim("~/Desktop/Landrace_Analysis/Worked_Datasets/GeneticMap_iSelect_9k.txt")

loc.df <- `803_landraces_KML`
geno.df <- Land_6152_SNPs_AB

options(stringsAsFactors = FALSE)

#replace character alleles, with integer genotypes
for(i in (1:ncol(geno.df))){
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

#### E or W 
    # West is 1, East is 2
  loc.df$side <- ifelse(test = loc.df$Longitude <48, yes = "1", no = "2")
    # making a temp df to store the names and the side column
  locside.df <- loc.df[,c("Accession.ID","side")]
     # merging onto the genotypes 
  locgenoEW.df <- merge(x = locside.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
     #basic.stats
  EW_basicstats <- basic.stats(locgenoEW.df[, -1], diploid = FALSE)
  write.csv(EW_basicstats$perloc,file="EW_basicstats.csv")
     #graph outlier SNPs
  EW_basicstats <- read.csv("~/Desktop/EW_basicstats.csv")
  EW_Fst <- EW_basicstats[,c("X","Fst")]
  SNP_by_chro <- GeneticMap_iSelect_9k[,c("SNP","cm","Chromosome")]

  SNP_FST.df <- merge(x = EW_Fst, y = SNP_by_chro, by.x = "X", by.y = "SNP")

  with(SNP_FST.df, plot(x = jitter(Chromosome,2),y = Fst,
                        cex = .75, col = ifelse(Fst > .5, "red", "black"), pch = ifelse(Fst > .5, 19, 21),
                        xlab = "Chromosome", ylab = "Fst value", main = "Outlier Fst Values"))
  abline(h=.5, lty = 6, lwd = .5, col = "red")
    #table of outliers
  Outliers <- subset(SNP_FST.df, Fst > .5)
  tgeno <- read.csv("~/Desktop/Landrace_Analysis/Worked_Datasets/tgeno.csv")
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
  
  #### Making a column for 2 v 6 row barley
  Rowtype.df <- barley_compare_merged[,c(2,11)]
  Rowtype2.df <- subset(Rowtype.df, SPIKEROW == 2, select = c(Accession.ID, SPIKEROW))
  Rowtype6.df <- subset(Rowtype.df, SPIKEROW == 6, select = c(Accession.ID, SPIKEROW))
  RowtypeC.df <- rbind(Rowtype2.df,Rowtype6.df)
  locRowtype.df <- merge(x = RowtypeC.df,y = geno.df, by.x = "Accession.ID", by.y = "X")
  locRowtype.df$Accession.ID <- 1:nrow(locRowtype.df)
  # 2:2; 6:1; other:NA
  locRowtype.df$SPIKEROW <- ifelse(test = locRowtype.df$SPIKEROW == "2", yes="2", no ="1" )
  row_basicstats <- basic.stats(locRowtype.df[, -1], diploid = FALSE)
  row_basicstats #FST: 0.0766
  #row_pp.fst <- pp.fst(dat=locRowtype.df,diploid = FALSE)
  
#### Making a column for spring vs winter
  SpringLR <- subset(mergedAccessions, HABIT == "S", select = c(Accession.ID,HABIT))
  WinterLR <- subset(mergedAccessions, HABIT == "W", select = c(Accession.ID,HABIT))
  HabitC.df <- rbind(SpringLR,WinterLR)
  locHabit.df <- merge(x = HabitC.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  locHabit.df$Accession.ID <- 1:nrow(locHabit.df)
    # Winter:1, Spring:2
  locHabit.df$HABIT <- ifelse(test = locHabit.df$HABIT == "W", yes="1", no="2")
  SpWi_basicstats <- basic.stats(locHabit.df[, -1], diploid = FALSE)
  SpWi_basicstats #FST = 0.0231
  

