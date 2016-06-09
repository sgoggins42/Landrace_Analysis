# No factors please
options(stringsAsFactors = FALSE)

# Import dataset
#   Genotype data
    geno.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/Fumi.Cleaned.Genotypes.csv")
#   Save Sample names and SNP names:
#   Samples:
    Samples.dat <- as.data.frame(geno.dat$taxa)
#   SNPs
    SNPs.dat <- as.data.frame(colnames(geno.dat))    
    SNPs.dat <- as.data.frame(SNPs.dat[-1,])
#   Genotype matrix only:
    geno.dat <- geno.dat[,-1]
# Turn genotype calls into minor allele counts:
#   Function to replace genotypes with A or B
    to.Assignment <- function(dat){
      dat[dat == 0] <- "A"
      dat[dat == 1] <- NA
      dat[dat == 2] <- "B"
      return(dat)
    }
#   Run function:
    geno.dat <- to.Assignment(dat = geno.dat)   
# Find columns where A is major allele:
  major <- names(which(apply(geno.dat, 2, function(x) length(which(x == "A")) / nrow(geno.dat)) >= 0.5))
# Find columns where A is minor allele:
  minor <- names(which(apply(geno.dat, 2, function(x) length(which(x == "A")) / nrow(geno.dat)) < 0.5))
# Subset these datasets
  major.dat <- geno.dat[,major]
  minor.dat <- geno.dat[,minor]
# Function from A major to 0, B minor to 1:
  to.A.major <- function(dat){
    dat[dat == "A"] <- as.numeric(0)
    dat[dat == "B"] <- as.numeric(1)
    return(dat)
  }
#  Run function
  major.dat <- to.A.major(dat = major.dat)     
# Function from A minor to 1, B major to 0:
  to.A.minor <- function(dat){
    dat[dat == "A"] <- as.numeric(1)
    dat[dat == "B"] <- as.numeric(0)
    return(dat)
  }
#  Run function
  minor.dat <- to.A.minor(dat = minor.dat)    
# Combine dataframes into one matrix of major allele counts:
  geno.dat <- as.data.frame(cbind(major.dat, minor.dat))
# Make data numeric:
  geno.dat <- as.data.frame(apply(geno.dat, 2, as.numeric))
  
# Plot the minor allele frequency:
    h <- hist(colSums(geno.dat, na.rm = TRUE)/nrow(geno.dat),
              breaks = 40)
    h$density <- h$counts/sum(h$counts)    
    plot(h, freq = FALSE,
         main = "MAF",
         xlab = "MAF",
         ylab = "Density")    
### What is going on here?


####
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
  

