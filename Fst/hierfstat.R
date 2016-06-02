#install.packages("hierfstat") #depends on 'gtools' and 'ade4'
library("hierfstat")

# No factors please
options(stringsAsFactors = FALSE)

# import datasets
#   Location data
  loc.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")
#   Genotype data
  geno.dat <- read.delim("~/Landrace_Analysis/Worked_Datasets/Land_6152_SNPs_AB.txt")
#   Mapping data
  GeneticMap_iSelect_9k <- read.delim("~/Landrace_Analysis/Worked_Datasets/GeneticMap_iSelect_9k.txt")

# Function to replace allles with integery genotypes
# AA = 11, and BB = 22, AB/BA = 12
  to.Assignment <- function(dat){
    dat[dat == "AA"] <- as.numeric(11)
    dat[dat == "BB"] <- as.numeric(22)
    dat[dat == "AB"] <- as.numeric(12)
    dat[dat == "BA"] <- as.numeric(12)
  return(dat)
  }
# Run function:
  geno.dat<- to.Assignment(dat = geno.dat)

# Make values numeric type, and add back in Sample names and rownames
#   Samples are the accession names
  Samples <- geno.dat[,1]
#   Samples are removed when making data numeric
  geno.dat <- as.data.frame(apply(geno.dat[,-1], 2, as.numeric))
#   Add Samples back in
    rownames(geno.dat) <- Samples
# Filtering:
# Remove mon-allelic SNPs
#   The column sum should be greater than if it only had 11 and less than if it only had 22
  geno.dat <- geno.dat[ , colSums(geno.dat, na.rm = TRUE) > (11 * length(geno.dat[,1])) & colSums(geno.dat, na.rm = TRUE) < (22 * length(geno.dat[,1]))]
# Remove SNPs with > 10% NA's
# Names of those columns that have 90% or greater data
  Data.NA <- names(which(apply(geno.dat, 2, function(x) sum(length(which(is.na(x))))) <= (0.1 * length(geno.dat[,1]))))
# Subset only these SNPs
  geno.dat <- geno.dat[ , Data.NA]
# Omit samples that have > 10% heterozygosity
  Data.hets <- names(which(apply(geno.dat, 1, function(x) sum(length(which(x == 12)))) <= (0.01 * length(geno.dat[,1]))))
# subset only these Samples
  geno.dat <- geno.dat[Data.hets, ]

####
# Trial:
  # scratch <- geno.dat
# Bring back the AA's and BB's
  # for(i in (1:ncol(scratch))){
  #   if((length(which(scratch[,i] == '11'))) > (length(which(scratch[,i] == '22')))){
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '11', "AA")
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '22', "BB")
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '12', "AB")
  # }
  # else{
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '11', "BB")
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '22', "AA")
  #   scratch[,i] <- replace(scratch[,i], scratch[,i] == '12', "AB")
  # }
  # }
# Make genotypes
# package  
    # library(genetics)
# command
    # scratch <- makeGenotypes(scratch, sep = '')
# color pallet
    #rgb.palette <- colorRampPalette(rev(c("blue","orange","red")),space = "rgb")
# Linkage map
    #My.map <- LDheatmap(scratch, LDmeasure = "r", distances = "genetic", title = "LD Heat Map", color = rgb.palette(18), flip=TRUE)
  
# Allele frequencies:
#   Names of SNPs with 11 as the minor allele  
    minor <- names(which(apply(geno.dat, 2, function(x) length(which(x == 11))) / nrow(geno.dat) <= 0.5))
#   Total frequncies of alleles    
    freqs <- apply(geno.dat, 2, function(x) length(which(x == 11))) / nrow(geno.dat) 
#   Grab the minor alleles
    freqs.plot <- freqs[minor]
#   Sort them by AF
    freqs <- sort(freqs.plot, decreasing = TRUE)
#   Barplot
    # barplot(height = freqs)

####  
# Fst over sliding latitudinal gradient
  library(genetics)
  library(hierfstat)
# Import data:
#   Fumi's cleaned Genotype data
    # Fumi.Cleaned.Genotypes <- read.csv("~/Landrace_Analysis/Worked_Datasets/Fumi.Cleaned.Genotypes.csv")
#   Re-location data
    # loc.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")
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
  grouped.genos.df <- merge(x = tmp.df, y= Genos, by.x = "Accession.ID", by.y = "taxa")
  grouped.genos.df$Lat <- apply(X = grouped.genos.df, MARGIN = 1, FUN = as.integer)
    # basic.stats, removing names from the calculations
  These.Stats <- basic.stats(grouped.genos.df[, -1], diploid = FALSE)
  These.F.Stats <- These.Stats$perloc
    # Report mean and range of top 5%
  This.mean <- data.frame(Mean = mean(These.F.Stats$Fst, na.rm = TRUE))
  This.overall <- data.frame(Overall = These.Stats$overall)
  This.Split <- data.frame(Lat = this.lat)
  Vals <- cbind(This.mean, This.overall, This.Split)
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

geno.dat <- Fumi.Cleaned.Genotypes

# Print graphs and stats
  # pdf(file = "Latitude Graphs")
  Fst.by.Lat <- sapply(X = Lats, FUN = to.Lats, dat = loc.dat, Genos = geno.dat)
  # dev.off()
  Fst.by.Lat

  # Get top Fst vals
  #SNP.Names <- data.frame(X = rownames(These.F.Stats))
  #FST.vals <- data.frame(Fst = These.F.Stats$Fst)
  #Top.dat.1 <- cbind(SNP.Names, FST.vals)
  #Top.dat.2 <- subset(x = Top.dat.1, subset = Fst >= Top.Fst)

  #### E or W 
  # West is 1, East is 2
  loc.df$side <- ifelse(test = loc.df$Longitude <48, yes = "1", no = "2")
  # making a temp df to store the names and the side column
  locside.df <- loc.df[,c("Accession.ID","side")]
  # merging onto the genotypes 
  locgenoEW.df <- merge(x = locside.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  #basic.stats
  EW_basicstats <- basic.stats(locgenoEW.df[, -1], diploid = FALSE)
  write.csv(EW_basicstats$perloc,file="EW_basicstats.csv")
  #graph outlier SNPs
  
  #### Making a column for 2 v 6 row barley
  #Rowtype.df <- barley_compare_merged[,c(2,11)]
  #Rowtype2.df <- subset(Rowtype.df, 
  #                      SPIKEROW == 2, 
  #                      select = c(Accession.ID, SPIKEROW))
  #Rowtype6.df <- subset(Rowtype.df, 
  #                        SPIKEROW == 6, 
  #                        select = c(Accession.ID, SPIKEROW))
  #RowtypeC.df <- rbind(Rowtype2.df,Rowtype6.df)
  #locRowtype.df <- merge(x = RowtypeC.df, y = geno.df, 
  #                       by.x = "Accession.ID", by.y = "X")
  #locRowtype.df$Accession.ID <- 1:nrow(locRowtype.df)
  # 2:2; 6:1; other:NA
  #locRowtype.df$SPIKEROW <- ifelse(test = locRowtype.df$SPIKEROW == "2", yes="2", no ="1" )
  #row_basicstats <- basic.stats(locRowtype.df[, -1], diploid = FALSE)
  #row_basicstats #FST: 0.0766
  
  #### Making a column for spring vs winter
  #SpringLR <- subset(mergedAccessions, HABIT == "S", select = c(Accession.ID,HABIT))
  #WinterLR <- subset(mergedAccessions, HABIT == "W", select = c(Accession.ID,HABIT))
  #HabitC.df <- rbind(SpringLR,WinterLR)
  #locHabit.df <- merge(x = HabitC.df, y=geno.df, by.x = "Accession.ID", by.y = "X")
  #locHabit.df$Accession.ID <- 1:nrow(locHabit.df)
  # Winter:1, Spring:2
  #locHabit.df$HABIT <- ifelse(test = locHabit.df$HABIT == "W", yes="1", no="2")
  #SpWi_basicstats <- basic.stats(locHabit.df[, -1], diploid = FALSE)
  #SpWi_basicstats #FST = 0.0231
  

