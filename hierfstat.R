install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

#import genotype dataset "Land_6152_SNPs_AB"

#basic.stats(data,diploid=TRUE,digits=4) #Fst
?basic.stats
col1 <- as.integer(Land_6152_SNPs_AB[,1])
nrow(Land_6152_SNPs_AB) #803
ncol(Land_6152_SNPs_AB) #6142
class(Land_6152_SNPs_AB[3,3]) #currently working with factors, need integers
#AA: 1
#BB: 2
#AB: 3
#is.na("NA)

#Success 
for(i in 1:ncol(Land_6152_SNPs_AB)){
  {
    if(i=="AA"){
      i <- 1
    }
    if(i=="BB"){
      i <- 2
    }
    if(i=="AB" & i=="BA"){
      i <- 3
    }
    is.na("NA") 
  }
}

View(Land_6152_SNPs_AB)

col1 <- Land_6152_SNPs_AB[1,]
basic.stats(col1, diploid=FALSE) #Error in 1 - sp2 - sHo/2/n : non-conformable arrays or Error in n/(n - 1) * (1 - sp2) : non-conformable arrays

basic.stats(Land_6152_SNPs_AB[,1],diploid=FALSE) #haploid?
?basic.stats
