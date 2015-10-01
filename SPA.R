#Select the landraces only
for(i in (1:ncol(geno.df))){
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AA",1)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BB",2)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="AB",NA)
  geno.df[,i]<-replace(geno.df[,i],geno.df[,i]=="BA",NA)
}

latlong_land<-loc.df[(loc.df[,1] %in% geno.df[,1]),]
dim(latlong_land)
head(latlong_land)
ZERO<-rep(0,803)
lacation_file<-cbind(as.data.frame(latlong_land[,1]),as.data.frame(latlong_land[,2]),as.data.frame(latlong_land[,3]),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(ZERO))
head(lacation_file)
colnames(lacation_file)[1] <- "Sample.ID"
colnames(lacation_file)[2] <- "Lat"
colnames(lacation_file)[3] <- "Lon"

minorAllele<- function(dat){
  count_A<-length(which(dat == '1'))
  count_B<-length(which(dat == '2'))
  #if A is > than B then B is the minor allele
  minor_allele<-if (count_A > count_B) {'B'} else {'A'}
  return(minor_allele)
}

MINOR_allele<-as.data.frame(geno.df)

##major_allele<-NULL
for (i in 1:(dim(MINOR_allele)[1])) {
  major_allele[i]<-if(MINOR_allele[i,] == '1'){'B'} else {'A'}
}

MAJOR_allele<-as.data.frame(major_allele)
