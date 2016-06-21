# No factors please
options(stringsAsFactors = FALSE)

# Import data set:
Bioclim <- read.csv("~/Landrace_Analysis/Worked_Datasets/Bioclim.csv")

# Get PC of Lat/Long, altitude, and bioclime data
Bioclim.pca <- prcomp(Bioclim[,c(2:23)], scale = TRUE)

#Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
#   Eigen vector
    loadings <- Bioclim.pca$rotation
#   Sd
    sdev <- Bioclim.pca$sdev
#   Eigen value
    varvector <- Bioclim.pca$sdev^2
#   Arrows to graph: (loadings*sdev)
    var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))


# Plot the correlation circle
  a <- seq(0, 2*pi, length = 100)
  plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC1",  ylab = "PC2",
      main = "Correlation between environmental variables \n and principal components",
      bg = gray(c(0.8)))

# Add active variables
# var.coord[,1]: PC1
# var.coord[,2]: PC2
  arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2,
       col = c("red","grey","grey","grey","grey","grey","grey","grey","red","grey","grey","grey",
               "grey","red","red","grey","grey","grey","grey","grey","grey","grey"))

# Add labels
  text(var.coord+.05, labels=rownames(var.coord), cex = 0.55, adj=1)
# Add lines to graph
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)

#The proportion of the variance that each eigenvector represents can be calculated by dividing the eigenvalue corresponding to that eigenvector by the sum of all eigenvalues.
  # Latitude
  # Bio6: 'Min Temperature of Coldest Month'
  # Bio11: 'Mean Temperature of Coldest Quarter'
  # Bio12: 'Annual Precipitation'


# princomp 3D graph: shows no trends
  pc <- princomp(Bioclim[,-1], cor=TRUE, scores=TRUE)
  plot3d(pc$scores[,1:3], col = c("red","blue","black"))

# Graph of GWAS results:
  GWAS.dat <- Bioclim[,c("Taxa","bio19")]
  GWAS.to.merge <- geno.dat
  GWAS.dat <- merge(x = GWAS.dat, y = GWAS.to.merge, by.x = "Taxa", by.y = "Accession.ID")
  t.test(GWAS.dat$bio19~GWAS.dat$SCRI_RS_219749)
  
  to.T.dist <- function(x){
    t.dat <- as.data.frame(t.test(GWAS.dat$bio19~x)$statistic))
    t.dat <- data.frame(t = t.dat$)
  }
  
  