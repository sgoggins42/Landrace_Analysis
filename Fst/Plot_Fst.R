# No factors please
options(stringsAsFactors = FALSE)

# Transform data by running hierfstat.R script
# source(hierfstat.R)
#   geno.dat: transformed for minor allele count
#   Samples.dat: SNP names from geno.dat
#   EW_perloc: output from hierfstat
#   Map.Locations: transformed from Genetic_iselect_map

# 1. Plot the minor allele frequency:
#   Major alleles only
    Major.alleles <- apply(geno.dat, 2, function(x) length(which(x == 0)))    
#   Minor alleles only  
    Minor.alleles <- apply(geno.dat, 2, function(x) length(which(x == 1)))    
#   Get frequencies
    Allele.freqs <- (Minor.alleles) / (Major.alleles + Minor.alleles)    
#   Plot
    pdf(file = "~/Landrace_Analysis/MAF.pdf")
    h <- hist(Allele.freqs,
              breaks = 40)
    h$density <- h$counts/sum(h$counts)    
    plot(h, freq = FALSE,
         main = "MAF",
         xlab = "MAF",
         ylab = "Density")    
    dev.off()
  
#2. Fst by Latitude analysis
#   Merge Fst with genetic mapping infor
    Map.Locations <- merge(x = EW_perloc, y = Map.Locations, by.x = "SNPs", by.y = "SNP")
#   Make chromosomes actually numeric:
    to.Chrom <- function(dat){
      dat[dat == "1H"] <- as.numeric(1)
      dat[dat == "2H"] <- as.numeric(2)
      dat[dat == "3H"] <- as.numeric(3)
      dat[dat == "4H"] <- as.numeric(4)
      dat[dat == "5H"] <- as.numeric(5)
      dat[dat == "6H"] <- as.numeric(6)
      dat[dat == "7H"] <- as.numeric(7)
      return(dat)
    }
#  Run function
    Map.Locations <- to.Chrom(dat = Map.Locations)   
#   Graph outlier SNPs
    library(qqman)
#   Change to fit manhattan
#   Real columns: SNP, Chro, cm, Fst
    Map.Locations <- Map.Locations[,c(1,3,4,2)]
    colnames(Map.Locations) <- c("SNP","CHR","BP","P")
    Map.Locations$CHR <- as.numeric(Map.Locations$CHR)
#   SNPs in pericentromereic regioin
    snpsOfInterest.2 <- subset(x = Map.Locations, 
                               subset = BP >= 67.1 & BP <= 68.1, select = c(SNP, CHR))
    snpsOfInterest.2 <- subset(x = snpsOfInterest.2, 
                               subset = CHR == 2, select = SNP)
    snpsOfInterest.3 <- subset(x = Map.Locations, 
                               subset = BP >= 65.3 & BP <= 67.9, select = c(SNP, CHR))
    snpsOfInterest.3 <- subset(x = snpsOfInterest.3, 
                               subset = CHR == 3, select = SNP)
    snpsOfInterest.4 <- subset(x = Map.Locations, 
                               subset =BP >= 56.1 & BP <= 56.3, select = c(SNP, CHR))
    snpsOfInterest.4 <- subset(x = snpsOfInterest.4, 
                               subset = CHR == 4, select = SNP)
    snpsOfInterest.5 <- subset(x = Map.Locations, 
                               subset = BP >= 42.2 & BP <= 42.5, select = c(SNP, CHR))
    snpsOfInterest.5 <- subset(x = snpsOfInterest.5, 
                               subset = CHR == 5, select = SNP)
    snpsOfInterest.6 <- subset(x = Map.Locations, 
                               subset = BP >= 59.3 & BP <= 60.6, select = c(SNP, CHR))
    snpsOfInterest.6 <- subset(x = snpsOfInterest.6, 
                               subset = CHR == 6, select = SNP)
    snpsOfInterest.7 <- subset(x = Map.Locations, 
                               subset = BP >= 80.3 & BP <= 81.8, select = c(SNP, CHR))
    snpsOfInterest.7 <- subset(x = snpsOfInterest.7, 
                               subset = CHR == 7, select = SNP)
    snpsOfInterest <- rbind(snpsOfInterest.2, 
                            snpsOfInterest.3, 
                            snpsOfInterest.4, 
                            snpsOfInterest.5,
                            snpsOfInterest.6, 
                            snpsOfInterest.7)
    snpsOfInterest <- as.array(snpsOfInterest$SNP)
#   Figure
    png(
      "/Users/dgrossen/Landrace_Analysis/Fst/Fst.by.Lat.png",
      width     = 3,
      height    = 2,
      units     = "in",
      res       = 1200,
      pointsize = 4
    )
   manhattan(x = Map.Locations, logp = FALSE,
                ylab = "FST",
                main = "FST by Latitude",
                highlight = snpsOfInterest)
   dev.off()
  
  
  