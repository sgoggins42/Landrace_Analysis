#SPA genotype input
geno.df <- Land_6152_SNPs_AB

for(i in (1:ncol(geno.df))){
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

# replace minor allele with 2, major allele with 0
for(i in (1:ncol(geno.df))){
  if((length(which(geno.df[,i] == '2'))) > (length(which(geno.df[,i] == '1')))){
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '2', 0)
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '1', 2)
  }
  else{
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '2', 2)
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '1', 0)
  }
}

spa_geno.df <- geno.df

#save populations names
Population_names <- spa_geno.df$X

#remove population names from dataset
spa_geno.df <- spa_geno.df[,-1]

#remove SNPs that arent biallelic
spa_geno.df <- spa_geno.df[, colSums(spa_geno.df == 0,na.rm=T) > 0  & colSums(spa_geno.df == 2,na.rm=T) > 0  ] 

#save SNP names
SNP_names <- colnames(spa_geno.df)
SNP_names <- as.data.frame(SNP_names)

write.table(spa_geno.df, sep = " ", file = "spa_geno_original", quote = FALSE, row.names = FALSE, col.names = FALSE)

