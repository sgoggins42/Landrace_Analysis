# Loc file output
    spa_loc <- read.delim("~/Desktop/Landrace_Analysis/SPA/spa_loc", stringsAsFactors=FALSE)
  #name columns
    cols <- c("Family.ID", "Individual.ID", "Paternal.ID", "Maternal.ID", "Sex", "Phenotype", "loc.dim1", "loc.dim2")
    names(spa_loc) <- cols

# SPA value output    
    spa_model <- read.delim("~/Desktop/Landrace_Analysis/SPA/spa_model", header=FALSE, stringsAsFactors=FALSE)
  #name columns
    cols2 <- c("Chromosome", "SNP.ID", "Genetic.distance","BP_Pos","Minor_Allele", "Major_Allele","Dim.1","Dim.2","Dim.3","SPA_Score")
    names(spa_model) <- cols2
    spa_model$SNP.ID <- 1:nrow(spa_model)
    
# Identify 95th percentile
    SPA_range <- as.data.frame(range(spa_model$SPA_Score) * .95)
    plot(x = spa_model$SNP.ID, y = spa_model$SPA_Score, main = "SPA values", xlab = "SNP", ylab = "SPA score", col = ifelse(test = spa_model$SPA_Score > SPA_range[2,], yes = "red", no = "blue"), log = "y")
  #make table out of outlier SPA values
    table(spa_model$SPA_Score > SPA_range[2,])
    high_SPA <- subset(spa_model, spa_model$SPA_Score > SPA_range[2,])

#Put back SNP names to number place holders from original spa input file
    #something about SNP.Names <- original.spa[,0]
    SNP.Names <- colnames(SNP.Names)
    SNP.Names<- as.data.frame(SNP.Names)
    high_SPA <- cbind(high_SPA, SNP.Names)
    # Reorganize: high_SPA <- catch[,c(10,2:9)]

#Identify SNP regions on chromosomes and such
    GeneticMap_iSelect_9k <- read.delim("~/Desktop/Landrace_Analysis/Worked_Datasets/GeneticMap_iSelect_9k.txt", stringsAsFactors=FALSE)
    GeneticMap <- GeneticMap_iSelect_9k
    GeneticMap <- GeneticMap[,c("SNP", "cm", "Chromosome")]
    outlier_SPA <- merge(x = high_SPA, y = GeneticMap, by.x = "SNP.Names", by.y = "SNP")
    colnames(outlier_SPA)[colnames(outlier_SPA)=="Chromosome.y"] <- "Chr"

#Add Fst info to dataframe
    outlier_SPA <- merge(x = outlier_SPA, y = EW_Fst, by.x = "SNP.Names", by.y = "X")
    #outlier_SPA <- outlier_SPA[,c(1:3,11,4:10)]

#Zhou's range on chromosome 2 is 67-74 cm
#out numbers are different, but both in range for one SNP.
Zhou_outlierSPA <- read.table("~/Desktop/Landrace_Analysis/SPA/Zhou_outlierSPA.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
Zhou_outlierSPA <- merge( x = Zhou_outlierSPA, y = outlier_SPA, by.x = "V1", by.y = "SNP.Names")
Zhou_outlierSPA$Z_cm <- 72.99
Zhou_outlierSPA <- Zhou_outlierSPA[, c(1:6, 12, 7:11)]


loc.df <- `803_landraces_KML`
geno.df$side <- 0
geno.spa <- geno.df
geno.spa$side <- ifelse(test = geno.spa$X11_10685 == 0, yes = "0", no = "2")
geno.spa <- geno.spa[,c("X","side")]
loc.spa <- merge(x = loc.df, y = geno.spa, by.x = "Accession.ID", by.y = "X")

library("maps")
data("wrld_simpl")
plot(wrld_simpl,xlim=c(-30,140),ylim=c(-8,64),col="lemonchiffon2",main="Distribution of Outlier Chromosome 2 SNP")
box()

with(loc.spa, points(x = Longitude,y = Latitude, cex = .75,
                      col = ifelse(side == "0", "red","blue")))
abline(v=48) #Zagaros 

legend("bottomleft", bg = "transparent",
       c("Zagros Mountains","Minor Allele", "Major Allele"),
       col=c("black","blue","red"), lwd=1, lty=c(1,NA,NA),
       pch = c(NA,21,21),cex = 0.75)