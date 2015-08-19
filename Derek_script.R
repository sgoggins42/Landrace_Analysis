#install.packages("heirfstat")
library("heirfstat")

setwd("~/Desktop/tmp/")
options(stringsAsFactors = FALSE)

loc.df <- read.csv("803_landraces_KML.csv")
geno.df <- read.csv("Land_6152_SNPs_AB.txt", sep ="\t", )

# reformat the genotype data to 1 and 2 

for(i in (1:ncol(geno.df))){
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

# making a column of E or W for location
 # West is 1, East is 2
loc.df$side <- ifelse(test = loc.df$Longitude <48, yes = "1", no = "2")

# making a temp df to store the names and the side column
locside.df <- loc.df[,c(1,6)]

# merging onto the genotypes 
locgeno.df <- merge(x = locside.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
head(locgeno.df[,1:3])
locgeno.df$Accession.ID <- 1:nrow(catch)
head(locgeno.df[,1:3])

EW_basicstats <- basic.stats(locgeno.df[, -1], diploid = FALSE)

