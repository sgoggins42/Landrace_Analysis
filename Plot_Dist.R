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
png(filename = "/Users/dgrossen/Landrace_Analysis/Landrace_plot.png")

this.map <- map("worldHires", xlim = c(-9,145), ylim=c(5,70), bg=gray(c(0.9)), fill=T, col= adjustcolor( "lemonchiffon", alpha.f = 0.8))

# Add points
  points(x = loc.dat$Longitude, y = loc.dat$Latitude, 
         pch = 18, cex = 0.75,
         col = ifelse(test = loc.dat$side ==1, yes = "red", no = "blue"),
         xlab = "Longitutde",
         ylab = "Latitude"
         )

  box(lwd=2)
  
  axis(side = c(1), at = seq(from = -5, to = 145, by = 5))
  axis(side = c(2), at = seq(from = 5, to = 70, by = 5))
  
  segments(x0 = 20, x1 = 72, y0 = 29, y1 = 29, lwd = 4, lty = 2, col = "chartreuse3")
  segments(x0 = 20, x1 = 72, y0 = 43, y1 = 43, lwd = 4, lty = 2, col = "chartreuse3")
  segments(x0 = 20, x1 = 20, y0 = 29, y1 = 43, lwd = 4, lty = 2, col = "chartreuse3")
  segments(x0 = 72, x1 = 72, y0 = 29, y1 = 43, lwd = 4, lty = 2, col = "chartreuse3")
  
dev.off()
  