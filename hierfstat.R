install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import genotype dataset "Land_6152_SNPs_AB"

#basic.stats(data,diploid=TRUE,digits=4) #Fst
#AA: 1
#BB: 2
#AB: NA

#import dataset
loc.df <- `803_landraces_KML`
geno.df <- Land_6152_SNPs_AB

#check data type
class(geno.df[3,4]) #needs to be character
options(stringsAsFactors = FALSE)

#replace character alleles, with integer genotypes
for(i in (1:ncol(geno.df))){
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
    geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

#### making a column of E or W for location
    # West is 1, East is 2
  loc.df$side <- ifelse(test = loc.df$Longitude <48, yes = "1", no = "2")
    # making a temp df to store the names and the side column
  locside.df <- loc.df[,c(1,6)]
  View(locside.df)
    # merging onto the genotypes 
  locgeno.df <- merge(x = locside.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  head(locgeno.df[,1:3])
  locgeno.df$Accession.ID <- 1:nrow(locgeno.df)
  head(locgeno.df[,1:3])
    #basic.stats
  EW_basicstats <- basic.stats(locgeno.df[, -1], diploid = FALSE)
  EW_basicstats #FST = 0.0802

#### Making a column for spring vs winter
  Habit.df <- barley_merged[ ,c(2,10)]
  head(Habit.df)
  locHabit.df <- merge(x = locHabit.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  head(locHabit.df[,1:3])
  locHabit.df$Accession.ID <- 1:nrow(locHabit.df)
    # Winter:1, Spring:2
  locHabit.df$HABIT <- ifelse(test = locHabit.df$HABIT == "W", yes="1", no="2")
  
  SpWi_basicstats <- basic.stats(locHabit.df[, -1], diploid = FALSE)
  SpWi_basicstats #FST = 0.0225
  
#Making a column for chr2 vs others
  genmap <- GeneticMap_iSelect_9k
  chr2<-subset(genmap,genmap$Chro == '2H')
  head(chr2)
  catchrows <- chr2[,c(1,2,7)]
  catchSNPs <- catchrows[,1]
  
  SNPs_head <- data.frame(x=catchSNPs)
  tgeno.df <- t(geno.df)
  View(tgeno.df)

  write.csv(tgeno.df,file="tgeno.csv")
  ?write.csv()
  
  SNPs_1.df <- merge(x = tgeno, y=SNPs_head, by.x = "X", by.y = "x")
  View(SNPs_1.df)
  
  
#Making a row for 2 v 5 v 6 row inflorescence 
  


  

