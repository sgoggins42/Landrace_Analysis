# No factors please
options(stringsAsFactors = FALSE)

# Import location data:
loc.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/803_landraces_KML.csv")

# Library - maps
library(maps)
library(mapdata)
library(mapplots)
library(scales)


loc.dat$side <- ifelse(test = loc.dat$Longitude <48, yes = "1", no = "2")
# Open plotting window
  #   quartz()
# Plot base map

#SCRI_RS_119778 #3

png(filename = "~/Landrace_Analysis/Landrace_plot.png")

this.map <- map("worldHires", xlim = c(-9,145), ylim=c(5,70), bg="transparent", fill=T, col= adjustcolor( "lemonchiffon", alpha.f = 0.8))

# Add points
  loc.scratch <- as.data.frame(merge(x = loc.dat, 
                                     y = geno.dat[,c("Accession.ID","X12_11190")],
                                     by = "Accession.ID") )
  
  points(x = loc.scratch$Longitude, y = loc.scratch$Latitude, 
         pch = 18, cex = 0.75,
         col = ifelse(test = loc.scratch$X12_11190 ==1, yes = "red", no = "blue"),
         xlab = "Longitutde",
         ylab = "Latitude"
  )

  box(lwd=2)
  
  axis(side = c(1), at = seq(from = -5, to = 145, by = 5))
  axis(side = c(2), at = seq(from = 5, to = 70, by = 5))

dev.off()
  