#!: Shawn Goggins - fst

#install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import dataset
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

#### Ana's populations
AllSNPs <- GeneticMap_iSelect_9k[,c("SNP","Chromosome")]
genoSNPs.df <- merge(x = tgeno, y=AllSNPs, by.x = "SNP.ID", by.y = "SNP") #chromosome at end

for(i in (1:ncol(Private_allele_ALLpops))){
  Private_allele_ALLpops[,i]<-replace(Private_allele_ALLpops[,i],Private_allele_ALLpops[,i]=="A",1)
  Private_allele_ALLpops[,i]<-replace(Private_allele_ALLpops[,i],Private_allele_ALLpops[,i]=="B",2)
}

#N_Mesopotamia private allele 
N_Mesopotamia.df <- subset(Private_allele_ALLpops, Population == "N_Mesopotamia", select = c("SNP", "Population"))
N_MesoSNPS.df <- merge(x=genoSNPs.df,y=N_Mesopotamia.df,by.x="SNP.ID",by.y="SNP")
s1 <- N_MesoSNPS.df$SNP.ID
N_MesoSNPS.df <- as.data.frame(t(N_MesoSNPS.df[,-1]))
colnames(N_MesoSNPS.df) <- s1
N_MesoSNPS.df <- N_MesoSNPS.df[1:803,]
N_MesoSNPS.df$Allele <- 1
N_MesoSNPS.df <- N_MesoSNPS.df[,c(ncol(N_MesoSNPS.df),1:(ncol(N_MesoSNPS.df)-1))]
#1
N_MesoSNPS.df$Allele <- 1:nrow(N_MesoSNPS.df)
#NM_basicstats <- basic.stats(N_MesoSNPS.df,diploid=FALSE)
#NM_basicstats
#write.csv(N_MesoSNPS.df,file = "N_MesoSNPS.df.csv")
#2
  #X12_30592 private allele 
    geno.df$Allele <- ifelse(test = N_MesoSNPS.df$X12_30592 == "2", yes = "1", no = "2")
    locNM.df <- geno.df[,c(1,ncol(geno.df),2:(ncol(geno.df)-1))]
    locNM.df$X <- 1:nrow(locNM.df)
    NMX12_30592_basicstats <- basic.stats(locNM.df[, -1], diploid = FALSE)
    NMX12_30592_basicstats #FST: 0.0388
    write.csv(NMX12_30592_basicstats$perloc, file = "NMX12_30592_perloc.csv")

S_Levant.df <- subset(Private_allele_ALLpops, Population == "S_Levant", select = c("SNP", "Population"))
S_LevantSNPS.df <- merge(x=genoSNPs.df,y=S_Levant.df,by.x="SNP.ID",by.y="SNP")

S_Desert.df <- subset(Private_allele_ALLpops, Population == "S_Desert", select = c("SNP", "Population"))
S_DesertSNPS.df <- merge(x=genoSNPs.df,y=S_Desert.df,by.x="SNP.ID",by.y="SNP")

N_Levant.df <- subset(Private_allele_ALLpops, Population == "N_Levant", select = c("SNP", "Population"))
N_LevantSNPS.df <- merge(x=genoSNPs.df,y=N_Levant.df,by.x="SNP.ID",by.y="SNP")

C_Asia.df <- subset(Private_allele_ALLpops, Population == "C_Asia", select = c("SNP", "Population"))
C_AsiaSNPS.df <- merge(x=genoSNPs.df,y=C_Asia.df,by.x="SNP.ID",by.y="SNP")


#### making a column of E or W for location
    # West is 1, East is 2
  loc.df$side <- ifelse(test = loc.df$Longitude <48, yes = "1", no = "2")
    # making a temp df to store the names and the side column
  locside.df <- loc.df[,c(1,6)]
    # merging onto the genotypes 
  locgenoEW.df <- merge(x = locside.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  locgenoEW.df$Accession.ID <- 1:nrow(locgenoEW.df)
    #basic.stats
  EW_basicstats <- basic.stats(locgenoEW.df[, -1], diploid = FALSE)
  EW_basicstats #FST = 0.0802
  #EW_ppfst <- pp.fst(dat = locgenoEW.df, diploid=FALSE)
  
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
  

 
