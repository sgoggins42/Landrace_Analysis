#install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import datasets
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

write.csv(geno.df, file = "geno.df")

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
  

