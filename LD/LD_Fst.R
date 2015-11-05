#run hierfstat.R to get: genodf
#source("hierfstat.R")
#run East_West_biostats.R to get: chro.2.inversion
#source("East_West_biostats.R")

library('LDheatmap')
library('genetics')

#Select SNPs that are outliers according to "chro.2.inversion"
chromosome2.outliers <- geno.df[,c("X", "SCRI_RS_185513", "SCRI_RS_232669", "X11_10317", "X11_11178", "X12_30561", "X12_30582")]

#only use samples with complete data
complete_chro2.df <- na.omit(chromosome2.outliers)

#find minor and major alleles
low_freq <- function(dat){
  A <- dat [dat == 1]
  B <- dat [ dat == 2]
  smaller <- min(c(A,B))
  return(c(length(A),length(B)))
}
#find major/minor allele: B
freq_SCRI_RS_185513 <- low_freq(complete_chro2.df$SCRI_RS_185513)
freq_SCRI_RS_232669 <- low_freq(complete_chro2.df$SCRI_RS_232669)
freq_X11_10317 <- low_freq(complete_chro2.df$X11_10317)
freq_X11_11178 <- low_freq(complete_chro2.df$X11_11178)
freq_X12_30561 <- low_freq(complete_chro2.df$X12_30561)
freq_X12_30582 <- low_freq(complete_chro2.df$X12_30582)

#Make dataframe LDheatmap friendly and run
genos <- complete_chro2.df[,-1]

#all SNPs are 1 = major, 2 = minor
for(i in (1:ncol(genos))){
  genos[,i]<-replace(genos[,i],genos[,i]=="1","AA")
  genos[,i]<-replace(genos[,i],genos[,i]=="2","BB")
}

genos <- makeGenotypes(genos, sep = '')

rgb.palette <- colorRampPalette(rev(c("blue","orange","red")),space = "rgb")

chro2.linkage <- LDheatmap(genos, LDmeasure = "r", distances = "genetic", title = "Chromosome 2 Fst Outliers", color = rgb.palette(18), flip=TRUE) #100% linkage
