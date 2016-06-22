# No factors
  options(stringsAsFactors = FALSE)
# Read in data
#   source(Geno.to.Fst.R) # Get EW_perloc or Fst results and Map.Locations for mapping SNPs
#   Read in other data or source(Plot_SPA.R) # To read in SPA results 
    spa_geno <- read.delim("~/Landrace_Analysis/SPA/Run_SPA/Fumis.Data/spa_geno_model", 
                           header=FALSE)
    SNP.names <- as.data.frame(colnames(geno.dat))[-1,]
    SNP.names <- as.data.frame(SNP.names)
    colnames(SNP.names) <- "SNPs"
    spa_geno$V2 <- SNP.names$SNPs
# Get top 1% 
#   Merge all
    Test.stats <- as.data.frame(merge(x = spa_geno[,c("V2","V10")], y = EW_perloc, 
                                    by.x = "V2", by.y = "SNPs"))
    colnames(Test.stats) <- c("SNP","SPA","FST")
#   Top Fst

