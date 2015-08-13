install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import genotype dataset "Land_6152_SNPs_AB"

#basic.stats(data,diploid=TRUE,digits=4) #Fst
#AA: 11
#BB: 22
#AB: 12

#import dataset
Land_6152_SNPs_AB <- read.delim("~/Desktop/FST/BarleyLandraces-master-2/Datasets/Land_6152_SNPs_AB.txt", stringsAsFactors=FALSE)
View(Land_6152_SNPs_AB)

#check data type
class(Land_6152_SNPs_AB[3,4]) #needs to be character

#replace character alleles, with integer genotypes
for(i in (1:ncol(Land_6152_SNPs_AB))){
    Land_6152_SNPs_AB[,i]<-replace(Land_6152_SNPs_AB[,i],Land_6152_SNPs_AB[,i]=="AA",1)
    Land_6152_SNPs_AB[,i]<-replace(Land_6152_SNPs_AB[,i],Land_6152_SNPs_AB[,i]=="BB",2)
    Land_6152_SNPs_AB[,i]<-replace(Land_6152_SNPs_AB[,i],Land_6152_SNPs_AB[,i]=="AB",NA)
    Land_6152_SNPs_AB[,i]<-replace(Land_6152_SNPs_AB[,i],Land_6152_SNPs_AB[,i]=="BA",NA)
}

View(Land_6152_SNPs_AB)

#add row
Land_6152_SNPs_AB$Patch <-1:nrow(Land_6152_SNPs_AB)

#shift new row to beginning
Land_genotypes <- Land_6152_SNPs_AB[,c(6154,1:6153)]
View(Land_genotypes)

#aligns East and West
for(i in Land_genotypes$X){
  Land_genotypes$X<-replace(Land_genotypes$X,Land_genotypes[i,2]==row.names(EastLR),1)
  Land_genotypes$X<-replace(Land_genotypes$X,Land_genotypes[i,2]==row.names(WestLR),0)
}

#run basic.stats on East v West v Other dataset
EvWvO <- basic.stats(Land_genotypes[,-1],diploid=FALSE)