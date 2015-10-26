data_loc <- `803_landraces_KML`
data_loc <- data_loc[,c(1:4)]
data_loc$v1 <- 1:nrow(data_loc)
data_loc$v3 <- 0
data_loc$v4 <- 0
data_loc$v5 <- 0
data_loc$v6 <- 0
data_loc<- data_loc[ ,c("Accession.ID","v1","v3", "v4", "v5", "v6", "Longitude","Latitude")]
data_loc$v1 <- ifelse(test = data_loc$Longitude > 48, yes = 1, no = 0)
data_loc$Accession.ID <- 1:nrow(data_loc)
write.table(x = data_loc, file = "data_loc", quote = FALSE, row.names = FALSE, col.names = FALSE)
?map

map2 <- get_googlemap(center = c(lon = 51, lat = 35), 
                      zoom = 3, scale = 2,
                      maptype = "hybrid", language = "en-EN", 
                      sensor = FALSE, messaging = FALSE, urlonly = FALSE,
                      filename = "ggmap_google_spa", color = "color")
ggmap(map2) +
  geom_point(aes(x = Longitude, y = Latitude), data = `803_landraces_KML`, color = "#CC0000")

install.packages("scatterplot3d")
with(outlier_SPA, {
  scatterplot3d(Dim.2,
                Dim.1,
                Dim.3,
                main = "SPA with Infered Location")
})




