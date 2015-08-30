#!/bin/Rscript

#08/29/15
#Written by Shawn Goggins for NSF-GRF proposal
#Graph to depict barley landraces in question of environmental associations
library("maps")
library("maptools")

mergedAccessions <- read.csv("~/Desktop/for_Shawn/mergedAccessions.csv", stringsAsFactors=FALSE)

quartz()
data("wrld_simpl")
plot(wrld_simpl,xlim=c(-30,140),ylim=c(-8,64),col="lemonchiffon2",main="Barley Landraces")
box()

with(mergedAccessions, points(x = Longitude,y = Latitude, cex = .75,
                              col = ifelse(HABIT == "S", "red","blue"), 
                              pch = ifelse( SPIKEROW == "2",21,22 )))
abline(v=48) #Zagaros 

legend("bottomleft", bg = "transparent",
       c("Zagaros Mountains","Spring (6) Lines", "Spring (2) Lines", "Winter (6) Lines", "Winter (2) Lines"),
       col=c("black","red","red","blue","blue"), lwd=1, lty=c(1,NA,NA,NA,NA),
       pch = c(NA,22,21,22,21),cex = 0.8)
