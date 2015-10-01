#   Test function
start <- function() {
  print("Welcome to makeKML.R")
  #inputData <- readline("Please specify the CSV file that you would like to extract inforamtion from:")
  print("Please specify the CSV file that you would like to extract inforamtion from:")
  inputData <- readLines(con = "stdin", 1)
  return(inputData)
}

#   Read the file
readFile <- function(inputData) {
  sampleInfo <- read.csv(file = inputData, header = TRUE)
  return(sampleInfo)
}

#   A function to install required packages
pkgTest <- function(package) {
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}

#   Package dependency list
pkgList <- c("XML", "maps")

#   Load the packages
batchInstall <- function(pkgList) {
  options(repos = c(CRAN = "http://cran.rstudio.com"))
  for(dep in pkgList) {
    pkgTest(dep)
  }
  lapply(X = pkgList, FUN = library, character.only = TRUE)
}

#   Define a function to make the KML file
makeKML <- function(lat, long, alt, name) {
  #   Start KML document
  kmlConnection <- xmlOutputBuffer(nameSpace = "", nsURI = "http://opengis.net/kml/2.2")
  kmlConnection$addTag("kml", close = FALSE) # kml tag
  kmlConnection$addTag("Document", close = FALSE) # Document tag
  for(i in 1:length(lat)) {
    #   Iterate over all places defined previously
    kmlConnection$addTag("Placemark", close = FALSE) # Placemark tag
    kmlConnection$addTag("name", name[i]) # Name of place
    kmlConnection$addTag("description", name[i])
    kmlConnection$addTag("Point", close = FALSE) # It's a point
    kmlConnection$addTag("coordinates", paste(long[i], lat[i], alt[i], sep = ",")) # The coordinates themselves
    kmlConnection$closeTag() # Point
    kmlConnection$closeTag() # Placemark
  }
  kmlConnection$closeTag() # Document
  kmlConnection$closeTag() # kml tag
  kmlConnection$value() # Return the values, spitting out kml to stdout
}

#   Extract information from the CSV file
pullColumns <- function(sampleInfo) {
  names <- sampleInfo[, 1]
  lattitudes <- sampleInfo[, 2]
  longitudes <- sampleInfo[, 3]
  altitudes <- sampleInfo[, 4]
  necessary <- list(names = names, lattitudes = lattitudes, longitudes = longitudes, altitudes = altitudes)
  return(necessary)
}

#   Do the work here
main <- function(sampleInfo = NULL, name = "") {
  # print("Welcome to makeKML.R")
  # #inputData <- readline("Please specify the CSV file that you would like to extract inforamtion from:")
  # print("Please specify the CSV file that you would like to extract inforamtion from:")
  # inputData <- readLines(con = "stdin", 1)
  if(is.null(sampleInfo)){ 
    inputData <- start()
    sampleInfo <- readFile(inputData = inputData)
  }
  batchInstall(pkgList = pkgList)
  necessary <- pullColumns(sampleInfo = sampleInfo)
  cat(makeKML(lat = necessary$lattitudes, long = necessary$longitudes, alt = necessary$altitudes, name = necessary$names), 
      file = paste(getwd(),"/coordinates",name,".kml",sep=""))  # takes stdout from kml and redirects (cat) to a file
}


#EastLR <- subset(x = `803_landraces_KML`, `803_landraces_KML`$Longitude > 48)
#main(EastLR, "_east")

#WestLR <- subset(x = `803_landraces_KML`, `803_landraces_KML`$Longitude < 48)
#main(WestLR, "_west")

main(<file>) #file should be in local directory as coordinates.kml
