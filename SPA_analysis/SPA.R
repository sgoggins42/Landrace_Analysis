#SPA genotype input
geno.df <- Land_6152_SNPs_AB

for(i in (1:ncol(geno.df))){
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

# replace minor allele with 1, major allele with 0
for(i in (1:ncol(geno.df))){
  if((length(which(geno.df[,i] == '2'))) > (length(which(geno.df[,i] == '1')))){
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '2', 0)
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '1', 1)
  }
  else{
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '2', 1)
    geno.df[,i] <- replace(geno.df[,i], geno.df[,i] == '1', 0)
  }
}

spa_geno.df <- geno.df

getwd()
write.csv(spa_geno.df, file = "Spa_input.df")
