# SPA value output    
spa_model <- read.delim("~/Desktop/Landrace_Analysis/SPA/Run_SPA/spa_model", header=FALSE, stringsAsFactors=FALSE)
  #name columns
    cols2 <- c("Chromosome", "SNP.ID", "Genetic.distance","BP_Pos","Minor_Allele", "Major_Allele","Dim.1","Dim.2","Dim.3","SPA_Score")
    names(spa_model) <- cols2
    spa_model$SNP.ID <- 1:nrow(spa_model)

#Separate SNPs by chromosome
SNP_by_chro <- GeneticMap_iSelect_9k[,c("SNP","cm","Chromosome")]

#Get SNP_names from write_spa_geno.R
SNP_by_SPA <- spa_model[,c("SNP.ID","SPA_Score")]
SNP_by_SPA$SNP.ID <- SNP_names$SNP_names

#Merge by SNPs
SNP_by_SPA <- merge(x = SNP_by_SPA, y = SNP_by_chro, by.x = "SNP.ID", by.y = "SNP")

# Identify 95th percentile
    SPA_range <- as.data.frame(range(spa_model$SPA_Score) * .95)
    Outlier_SPA <- subset(SNP_by_SPA, SPA_Score > SPA_range[2,1])
    
#Plot distribution of SPA scores
    with(SNP_by_SPA, plot(x = jitter(Chromosome,2),y = SPA_Score,
                          cex = .75, col = ifelse(SPA_Score > SPA_range[2,1], "red", "black"),
                          xlab = "Chromosome", ylab = "SPA Score", main = "Outlier SPA Values"))
    abline(h = SPA_range[2,1], lty = 6, lwd = .5, col = "red")
    