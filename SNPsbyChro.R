######### Makes a column for SNPs by chromosome
AllSNPs <- GeneticMap_iSelect_9k[,c("SNP","Chromosome")]
genoSNPs.df <- merge(x = tgeno, y=AllSNPs, by.x = "SNP.ID", by.y = "SNP") #chromosome at end

genoSNPs.df <- genoSNPs.df[,c(1,805,2:804)] #write.csv(genoSNPs.df,file="SNPs_by_chrom.csv")
View(genoSNPs.df[1:5,1:5])

chr1_geno <- subset(genoSNPs.df,Chromosome=="1")
n1 <- chr1_geno$SNP.ID
chr1_geno <- as.data.frame(t(chr1_geno[,-1]))
colnames(chr1_geno) <- n1
chr1_geno <- chr1_geno[-1,]
chr1_geno$chr <- 1
chr1_geno <-   chr1_geno[,c(ncol(chr1_geno),1:(ncol(chr1_geno)-1))] 
chr1_geno$ID <- row.names(chr1_geno)
chr1_geno <-   chr1_geno[,c(ncol(chr1_geno),1:(ncol(chr1_geno)-1))] 
View(chr1_geno[1:5,1:5])
write.csv(chr1_geno,file="chr1_geno")

chr2_geno <- subset(genoSNPs.df,Chromosome=="2")
n2 <- chr2_geno$SNP.ID
chr2_geno <- as.data.frame(t(chr2_geno[,-1]))
colnames(chr2_geno) <- n2
chr2_geno <- chr2_geno[-1,]
chr2_geno$chr <- 2
chr2_geno <-   chr2_geno[,c(ncol(chr2_geno),2:(ncol(chr2_geno)-1))] 
chr2_geno$ID <- row.names(chr2_geno)
chr2_geno <-   chr2_geno[,c(ncol(chr2_geno),1:(ncol(chr2_geno)-1))] 
View(chr2_geno[1:5,1:5])
write.csv(chr2_geno,file="chr2_geno")

chr3_geno <- subset(genoSNPs.df,Chromosome=="3")
n3 <- chr3_geno$SNP.ID
chr3_geno <- as.data.frame(t(chr3_geno[,-1]))
colnames(chr3_geno) <- n3
chr3_geno <- chr3_geno[-1,]
chr3_geno$chr <- 3
chr3_geno <-   chr3_geno[,c(ncol(chr3_geno),2:(ncol(chr3_geno)-1))] 
chr3_geno$ID <- row.names(chr3_geno)
chr3_geno <-   chr3_geno[,c(ncol(chr3_geno),1:(ncol(chr3_geno)-1))] 
View(chr3_geno[1:5,1:5])
write.csv(chr3_geno,file="chr3_geno")

chr4_geno <- subset(genoSNPs.df,Chromosome=="4")
n4 <- chr4_geno$SNP.ID
chr4_geno <- as.data.frame(t(chr4_geno[,-1]))
colnames(chr4_geno) <- n4
chr4_geno <- chr4_geno[-1,]
chr4_geno$chr <- 4
chr4_geno <-   chr4_geno[,c(ncol(chr4_geno),2:(ncol(chr4_geno)-1))] 
chr4_geno$ID <- row.names(chr4_geno)
chr4_geno <-   chr4_geno[,c(ncol(chr4_geno),1:(ncol(chr4_geno)-1))] 
View(chr4_geno[1:5,1:5])
write.csv(chr4_geno,file="chr4_geno")

chr5_geno <- subset(genoSNPs.df,Chromosome=="5")
n5 <- chr5_geno$SNP.ID
chr5_geno <- as.data.frame(t(chr5_geno[,-1]))
colnames(chr5_geno) <- n5
chr5_geno <- chr5_geno[-1,]
chr5_geno$chr <- 5
chr5_geno <-   chr5_geno[,c(ncol(chr5_geno),2:(ncol(chr5_geno)-1))] 
chr5_geno$ID <- row.names(chr5_geno)
chr5_geno <-   chr5_geno[,c(ncol(chr5_geno),1:(ncol(chr5_geno)-1))] 
View(chr5_geno[1:5,1:5])
write.csv(chr5_geno,file="chr5_geno")

chr6_geno <- subset(genoSNPs.df,Chromosome=="6")
n6 <- chr6_geno$SNP.ID
chr6_geno <- as.data.frame(t(chr6_geno[,-1]))
colnames(chr6_geno) <- n6
chr6_geno <- chr6_geno[-1,]
chr6_geno$chr <- 6
chr6_geno <-   chr6_geno[,c(ncol(chr6_geno),2:(ncol(chr6_geno)-1))] 
chr6_geno$ID <- row.names(chr6_geno)
chr6_geno <-   chr6_geno[,c(ncol(chr6_geno),1:(ncol(chr6_geno)-1))] 
View(chr6_geno[1:5,1:5])
write.csv(chr6_geno,file="chr6_geno")

chr7_geno <- subset(genoSNPs.df,Chromosome=="7")
n7 <- chr7_geno$SNP.ID
chr7_geno <- as.data.frame(t(chr7_geno[,-1]))
colnames(chr7_geno) <- n7
chr7_geno <- chr7_geno[-1,]
chr7_geno$chr <- 7
chr7_geno <-   chr7_geno[,c(ncol(chr7_geno),2:(ncol(chr7_geno)-1))] 
chr7_geno$ID <- row.names(chr7_geno)
chr7_geno <-   chr7_geno[,c(ncol(chr7_geno),1:(ncol(chr7_geno)-1))] 
View(chr7_geno[1:5,1:5])
write.csv(chr7_geno,file="chr7_geno")

library(plyr)
