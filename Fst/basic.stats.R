#install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

# No factors please
options(stringsAsFactors = FALSE)

# import datasets
#   Genotype data
  geno.dat <- read.delim("~/Landrace_Analysis/Worked_Datasets/Land_6152_SNPs_AB.txt")

# Make values numeric type, and add back in Sample names and rownames
#   Samples are the accession names
  Samples <- geno.dat[,1]
#   Samples are removed when making data numeric
  geno.dat <- as.data.frame(apply(geno.dat[,-1], 2, as.numeric))
#   Add Samples back in
  rownames(geno.dat) <- Samples
  
# Function to change integers to allele counts
  to.Assignment <- function(dat){
    dat[dat == 0] <- "BB"
    dat[dat == 2] <- "AA"
    return(dat)
  }
# Apply Function
  geno.dat <- to.Assignment(dat = geno.dat)
  
# Count AA's
minors <- apply(geno.dat, 2, function(x) (length(which(x == "AA"))))
# Which have AA's as the lesser count...are wrongly labled
oops.minors <- names(which( minors / nrow(geno.dat) < 0.5))
# Which have AA's as the greater count...are correctly labled
correct.minors <- names(which( minors / nrow(geno.dat) >= 0.5))
# Indexm data
oops.minors <- geno.dat[,oops.minors]
correct.minors <- geno.dat[,correct.minors]

# Write function to switch oops.minors
to.Opposite.Alleles <- function(dat){
  dat[dat == "AA"] <- as.integer(11)
  dat[dat == "BB"] <- as.integer(22)
  return(dat)
}
# Apply
oops.minors <- to.Opposite.Alleles(dat = oops.minors)

# Write function to change alleles to numbers for correct.minors
to.Alleles <- function(dat){
  dat[dat == "AA"] <- as.integer(22)
  dat[dat == "BB"] <- as.integer(11)
  return(dat)
}
# Apply
correct.minors <- to.Alleles(dat = correct.minors)

# combine datasets:
geno.for.fst <- cbind(oops.minors, correct.minors)
geno.for.fst <- as.data.frame(apply(X = geno.for.fst, MARGIN = 2, FUN = as.integer))
rownames(geno.for.fst) <- rownames(geno.dat)
geno.for.fst$X <- rownames(geno.for.fst)
####  
# Fst over sliding latitudinal gradient
# Packages
  library(genetics)
  library(hierfstat)
# Import data:
#   Location data
  loc.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")
# Data frame with a sliding window of latitutde values
  Lats <- as.list(seq(from = min(range(loc.dat$Latitude, na.rm = TRUE)), 
                                to = max(range(loc.dat$Latitude, na.rm = TRUE)),
                                by = 5))
# Function to get Fst: mean, overall (how are they different), latitude split, and top 1%SNPs:
to.Lats <- function(dat, Genos, this.lat) {  
    # dat = loc.dat
    # Genos = geno.dat
    # this.lat = Lats
    # Identify subpopulations based on latitude values above
  dat$Lat <- ifelse(test = dat$Latitude > this.lat, yes = 1, no = 2)
    # Making a temp df to store the names and the side column
  tmp.df <- dat[,c("Accession.ID","Lat")]
    # merging onto the genotypes 
  grouped.genos.df <- merge(x = tmp.df, y= Genos, by.x = "Accession.ID", by.y = "X")
  geno.names <- grouped.genos.df$Accession.ID
  grouped.genos.df$Lat <- apply(X = grouped.genos.df[,-1], MARGIN = 2, FUN = as.integer)
  rownames(grouped.genos.df) <- geno.names
    # basic.stats, removing names from the calculations
  These.Stats <- basic.stats(grouped.genos.df, diploid = FALSE)
  These.F.Stats <- These.Stats$perloc
    # Report mean and range of top 5%
  This.mean <- data.frame(Mean = mean(These.F.Stats$Fst, na.rm = TRUE))
  This.Split <- data.frame(Lat = this.lat)
  Vals <- cbind(This.mean, This.Split)
  # Graph
    # Get top 5%
  # Identify Fst Outliers by 95%
  #Grab <- 0.01 * length(These.F.Stats$Fst)
  #Sorted.Fst <- data.frame(sort(x = These.F.Stats$Fst, decreasing = TRUE))
  #colnames(Sorted.Fst) <- "Fst"
  #Top.Fst <- Sorted.Fst$Fst[Grab]
    # Histogram
  #hist(These.F.Stats$Fst,
  #     main = this.lat,
  #     xlab = "Fst",
  #     xlim = c(-.2, 1),
  #     breaks = 20,
  #     col = "grey50")
  #abline(v = Top.Fst, col = "red")
  # Return stats
    return(Stats = Vals)
} 


# Print graphs and stats
  # pdf(file = "Latitude Graphs")
  Fst.by.Lat <- sapply(X = Lats, FUN = to.Lats, dat = loc.dat, Genos = geno.for.fst)
  # dev.off()
  Fst.by.Lat

  # Get top Fst vals
  #SNP.Names <- data.frame(X = rownames(These.F.Stats))
  #FST.vals <- data.frame(Fst = These.F.Stats$Fst)
  #Top.dat.1 <- cbind(SNP.Names, FST.vals)
  #Top.dat.2 <- subset(x = Top.dat.1, subset = Fst >= Top.Fst)

  