# No factors Setting
options(stringsAsFactors = FALSE)
# Import data:
#   Genotype data
    geno.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/Fumi.Cleaned.Genotypes.csv")
#   Import location data:
    loc.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")
#   Import SNP location data:
    GeneticMap_iSelect_9k <- read.delim("~/Landrace_Analysis/Worked_Datasets/GeneticMap_iSelect_9k.txt", stringsAsFactors=FALSE)
# Import dataset
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
  

#2. Fst by Latitude analysis
#   Data frame with a sliding window of latitutde values
    Lats <- as.list(seq(from = min(range(loc.dat$Latitude, na.rm = TRUE)), 
                                to = max(range(loc.dat$Latitude, na.rm = TRUE)),
                                by = 5))
#   Column name change for Samples
    colnames(Samples.dat) <- "Accession.ID"
#   Add column to genotype matrix for merging later
    geno.dat <- as.data.frame(cbind(Samples.dat,geno.dat))
#   Import package
    library(hierfstat)
##  Sliding window Fst function
#   Function to get Fst: mean,, latitude split, and top 1% SNPs:
    to.Lats <- function(dat, Genos, this.lat) {  
    ## dat = loc.dat
    ## Genos = geno.dat
    ## this.lat = Lats
    ####
    #1. Perform Fst analysis:
    #   Identify subpopulations based on latitude values above
        dat$Lat <- ifelse(test = dat$Latitude > this.lat, yes = as.integer(1), no = as.integer(2))
    #   Make a temp df to store the names and the side column
        tmp.df <- dat[ ,c("Accession.ID","Lat")]
    #   Change column names for merging
        colnames(tmp.df) <- c("Accession.ID","Lat")
    #   Merge onto the genotypes 
        grouped.genos.df <- merge(x = tmp.df, y= Genos, by = "Accession.ID")
    #   hiefstat: basic.stats, removing names from the calculations
        These.Stats <- basic.stats(grouped.genos.df[, -1], diploid = FALSE)
    #   Save Fstats from perloc
        These.F.Stats <- These.Stats$perloc
    ####
    #2. Report mean and latiudinal split
    #   Mean:
        This.mean <- data.frame(Mean = mean(These.F.Stats$Fst, na.rm = TRUE))
    #   Split:
        This.Split <- data.frame(Lat = this.lat)
    #   Vals:
        Vals <- cbind(This.mean, This.Split)
    ####
    # 3. Graph the Fst distribution
    #   Means to add to graph
        graph.means <- as.list(This.mean$Mean)
    #   Histogram
        hist(These.F.Stats$Fst,
            main = "Distribution of Fst",
            xlab = "Fst",
            xlim = c(-.2, 1),
            breaks = 20,
            col = "grey50")
    #   Add means
        abline(v = graph.means, col = "red")
    #   Add legend
        legend("topright", legend = c("Latitude: ", this.lat, "Mean Fst:", graph.means),
             col = c(NA,NA,NA,"red"), lty = c(NA,NA,NA,1))
    ####
    # 4. Get top 1% SNPs per latitude group
       #Top.F <- data.frame(Fst = These.F.Stats$Fst, Accession.ID = SNPs.dat$SNPs.dat)
       #N <- length(this.lat)
       #vectorOfTables <- vector(mode = "list", length = N)
       #vectorOfTables[[this.lat]] <- Top.F[with(Top.F, 
        # order(Fst, na.last = TRUE, decreasing = TRUE)),]
    # Return stats
      return(Stats = Vals)
    } 
# End function ^
#   Run function and print graphs to output
    #pdf(file = "~/Landrace_Analysis/Latitude_Graphs.pdf")
    Fst.by.Lat <- sapply(X = Lats, FUN = to.Lats, dat = loc.dat, Genos = geno.dat)
    #dev.off()

#### E or W of Zagros Mts
#   West is 1, East is 2
    loc.dat$side <- ifelse(test = loc.dat$Longitude <48, yes = "1", no = "2")
    # making a temp df to store the names and the side column
    locside.df <- loc.dat[,c("Accession.ID","side")]
    # merging onto the genotypes 
    locgenoEW.df <- merge(x = locside.df, y=geno.dat, by = "Accession.ID")
    #basic.stats
    EW_basicstats <- basic.stats(locgenoEW.df[, -1], diploid = FALSE)
    EW_perloc <- EW_basicstats$perloc
    # Merge Fst with genetic map info
    #   Get SNP names from rownames to column
    EW_names <- as.data.frame(rownames(EW_perloc))
    colnames(EW_names) <- "SNPs"
    #   combine SNP names and FST
    EW_perloc <- as.data.frame(cbind(EW_names, EW_perloc$Fst))
    colnames(EW_perloc)[2] <- "Fst"
    #   Select map locations data
    Map.Locations <- GeneticMap_iSelect_9k[,c("SNP","Chro","cm")]
   