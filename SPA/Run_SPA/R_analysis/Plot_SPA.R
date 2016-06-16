options(stringsAsFactors = FALSE)

spa_geno <- read.delim("~/Landrace_Analysis/SPA/Run_SPA/Fumis.Data/spa_geno_model", header=FALSE)

spa_geno <- as.data.frame(spa_geno[,c("V2","V10")])

SNP.names <- as.data.frame(colnames(geno.dat))[-1,]
SNP.names <- as.data.frame(SNP.names)
colnames(SNP.names) <- "SNPs"

spa_geno$V2 <- SNP.names$SNPs


#   Select map locations data
Map.dat <- GeneticMap_iSelect_9k[,c("SNP","Chro","cm")]
#   Finally merge
Map.dat <- merge(x = spa_geno, y = Map.dat, by.x = "V2", by.y = "SNP")
# Make chromosomes actually numeric:
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
Map.dat <- to.Chrom(dat = Map.dat)   
#graph outlier SNPs
library(qqman)
# Change to fit manhattan
#   Real columns: SNP, Chro, cm, Fst
Map.dat <- Map.dat[,c(1,3,4,2)]
colnames(Map.dat) <- c("SNP","CHR","BP","P")
Map.dat$CHR <- as.numeric(Map.dat$CHR)


snpsOfInterest.2 <- subset(x = Map.dat, subset = BP >= 67.1 & BP <= 68.1, select = c(SNP, CHR))
snpsOfInterest.2 <- subset(x = snpsOfInterest.2, subset = CHR == 2, select = SNP)

snpsOfInterest.3 <- subset(x = Map.dat, subset = BP >= 65.3 & BP <= 67.9, select = c(SNP, CHR))
snpsOfInterest.3 <- subset(x = snpsOfInterest.3, subset = CHR == 3, select = SNP)

snpsOfInterest.4 <- subset(x = Map.dat, subset =BP >= 56.1 & BP <= 56.3, select = c(SNP, CHR))
snpsOfInterest.4 <- subset(x = snpsOfInterest.4, subset = CHR == 4, select = SNP)

snpsOfInterest.5 <- subset(x = Map.dat, subset = BP >= 42.2 & BP <= 42.5, select = c(SNP, CHR))
snpsOfInterest.5 <- subset(x = snpsOfInterest.5, subset = CHR == 5, select = SNP)

snpsOfInterest.6 <- subset(x = Map.dat, subset = BP >= 59.3 & BP <= 60.6, select = c(SNP, CHR))
snpsOfInterest.6 <- subset(x = snpsOfInterest.6, subset = CHR == 6, select = SNP)

snpsOfInterest.7 <- subset(x = Map.dat, subset = BP >= 80.3 & BP <= 81.8, select = c(SNP, CHR))
snpsOfInterest.7 <- subset(x = snpsOfInterest.7, subset = CHR == 7, select = SNP)

snpsOfInterest <- rbind(snpsOfInterest.2, snpsOfInterest.3, snpsOfInterest.4, snpsOfInterest.5,
                        snpsOfInterest.6, snpsOfInterest.7)

snpsOfInterest <- as.array(snpsOfInterest$SNP)

png(
  "/Users/dgrossen/Landrace_Analysis/SPA/SPA_Dist.png",
  width     = 3,
  height    = 2,
  units     = "in",
  res       = 1200,
  pointsize = 4
)

manhattan(x = Map.dat, logp = FALSE,
          ylab = "SPA",
          main = "Spatial Ancestry Analysis",
          highlight = snpsOfInterest)

dev.off()


