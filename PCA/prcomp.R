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
      main = "Correlation between variables and principal components")

# Add active variables
# var.coord[,1]: PC1
# var.coord[,2]: PC2
  arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2)

# Add labels
  text(var.coord, labels=rownames(var.coord), cex = 0.75, adj=1)
# Add lines to graph
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)

