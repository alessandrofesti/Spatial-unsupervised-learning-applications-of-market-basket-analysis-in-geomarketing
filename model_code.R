## Alessandro Festi, Spatial Association Rules for geomarketing purposes

# Import packages
library(geosphere)
library(shiny)
library(leaflet)
library(plyr) 
library(dplyr)
library(jsonlite)
library(arules)
library(MASS)
library(sp)
library(arulesViz)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(readxl)
library(reticulate)

use_python('/usr/bin/python3.8')

################ Python ################  

# Python   
# ```{python getting Latitude/longitude through the Mapbox APIs, echo=TRUE}
# 
#
# import os
# import pandas as pd
# import numpy as np 
# from mapbox import Geocoder
# import json
# 
# os.chdir('/home/fester/Scrivania')
# dataset = pd.read_csv("elenco_esercizi_commercio_in_sede_fissa_anno_2018.csv", sep = ';', header='infer', encoding='latin-1')
# dataset['quartiere_settore'] = dataset.ESERCIZIO_VIA+'  '+dataset.ESERCIZIO_CIVICO+' '+dataset.QUARTIERE+' Bologna'
# dataset['lat'] = float
# dataset['lon'] = float
# token = yourtoken
# geocoder = Geocoder(access_token=token)
# 
# def mapbox_geocode(dataset):
#   for i in range(len(dataset)):
#     try:
#       response = geocoder.forward(dataset.quartiere_settore[i])
#       ale = response.content
#       d = json.loads(ale)
#       coordinates = d["features"][0]['geometry']['coordinates']
#       dataset.iat[i,29] = coordinates[1]
#       dataset.iat[i,30] = coordinates[0]
#     except:
#       dataset.iat[i,29] = np.nan
#       dataset.iat[i,30] = np.nan
# return(dataset)
# 
# geocoded_dataset = mapbox_geocode(dataset)
# # write_csv2(geocoded, 'geocoded.csv')
# ```

########### End Python ############

################ R ################ 

setwd("...")
geocoded <- read_csv("geocoded.csv")
geocoded <- dplyr::distinct(geocoded, lat, lon, .keep_all = TRUE)

# Setting parameters
n_individuals <- 200
n_paths <- 50
n_positions <- 50
bologna_latitude <- 44.493674
bologna_longitude <- 11.342220
points <- data.matrix(cbind(rnorm(n_individuals)/200 + bologna_longitude, rnorm(n_individuals)/200 + bologna_latitude)) 
correlation <- 0.7

# Defining path-generator function
paths_gen <- function(points) {
  paths_one <- c()
  for (j in 1:nrow(points)) {  
    for (i in 1:n_paths) {           
      mu <- rep(0,n_positions) 
      Sigma <- matrix(correlation, nrow=n_positions, ncol=n_positions) + diag(n_positions)*.3
      rawvars <- mvrnorm(n=n_positions, mu=mu, Sigma=Sigma) 
      geo_points <- as.data.frame(cbind(rawvars[,i]/250 + as.double(points[j]), rawvars[,i]/250 + as.double(points[j+nrow(points)])))
      gg <- cbind(geo_points, i, j)
      paths_one <- rbind(paths_one, gg)
    }
  }  
  return(as.data.frame(paths_one))
}

path_vis <- as.data.frame(paths_gen(points))
path_vis <- plyr::rename(path_vis, c("V1"="lon", "V2"="lat", "i"="path", "j"="ind"))
head(path_vis,10)

path_v <- data.matrix(cbind(path_vis$lon,path_vis$lat))
places <- cbind(geocoded$lon, geocoded$lat)
places[is.na(places)] <- 0
min_dist <- 0.0005

# Defining function to assign individual positions to geocoded commercial activities
assign <- function(geo_points) {
  geo_pointsTOplaces <- list()
  for(j in 1:nrow(path_v)) {
    x <- spDistsN1(places,path_v[j,], longlat = FALSE)
    x <- as.data.frame(x)
    qq <- which(x == min(x), arr.ind = TRUE)
    qq <- as.data.frame(qq) 
    nearest <- sort(x[x>0],decreasing=F)[1]
    geo_pointsTOplaces$distance[j] <- nearest
    if (nearest < min_dist) { 
      geo_pointsTOplaces$lon[j] <- places[qq[1,1],1]
      geo_pointsTOplaces$lat[j] <- places[qq[1,1],2]
    } else {
      geo_pointsTOplaces$lon[j] <- NA
      geo_pointsTOplaces$lat[j] <- NA 
    }
  }
  return(geo_pointsTOplaces)
}

result <- as.data.frame(assign(geo_points))
result <- as.data.frame(cbind(result$lon, result$lat, result$distance, path_vis$path, path_vis$ind))
result <- plyr::rename(result, c("V1"="lon", "V2"="lat", "V3"="dist", "V4"="path", "V5"="ind"))
head(result)

# Delete non-matched observations
result <- na.omit(result)
head(result)

# Merging process
names(geocoded)[1] <- "IND"
index <- dplyr::select(geocoded,lat,lon, IND)
index <- plyr::rename(index, c("lat" = "lat", "lon" = "lon", "IND" = "ID"))
final <- merge(x=index, y=result, by.x = c("lon","lat"), by.y = c("lon","lat"))
final$itemset_id <- paste(final$ind,final$path)
head(final,10)

# Preparing for Association rules
transactionData <- ddply(final,c("itemset_id"),
                         function(final)paste(final$ID,

# write.csv2(transactionData,"transaction_data.csv", quote = FALSE, row.names = TRUE)
tr <- read.transactions('transaction_data.csv', format = 'basket', sep=',')
itemFrequencyPlot(tr,topN=20,type="relative",col=brewer.pal(8,'Pastel2'), main="Absolute Item Frequency Plot")

# Performing MBA through the apriori algorithm
association.rules <- arules::apriori(tr, parameter = list(supp=0.01, conf=0.3,maxlen=10))
inspect(association.rules[1:10])

# Making plots
plot(association.rules, method = "two-key plot")
plot(association.rules, method = "graph")
plot(association.rules[1:20], method = "graph")

# Working on categories
index_cat <- dplyr::select(geocoded,lat,lon, ATTIVITA_PREVALENTE_ESERCIZIO)
index_cat <- plyr::rename(index_cat, c("lat" = "lat", "lon" = "lon", "ATTIVITA_PREVALENTE_ESERCIZIO" = "Commercial type"))
final_cat <- merge(x=index_cat, y=result, by.x = c("lon","lat"), by.y = c("lon","lat"))
final_cat$itemset_id <- paste(final_cat$ind,final_cat$path)
final_cat <- na.omit(final_cat)
head(final_cat)

transactionData_cat <- ddply(final_cat,c("itemset_id"),
                             function(final_cat)paste(final_cat$`Commercial type`,
                                                      collapse = ","))

# write.csv(transactionData_cat,"transaction_data_cat.csv", quote = FALSE, row.names = TRUE)
tr_cat <- read.transactions('transaction_data_cat.csv', format = 'basket', sep=',')
itemFrequencyPlot(tr_cat,topN=10,type="relative",col=brewer.pal(15,'Pastel2'), main="Absolute Item Frequency Plot", horiz=TRUE)

association.rules_cat <- arules::apriori(tr_cat, parameter = list(supp=0.03, conf=0.5,maxlen=5))
inspect(association.rules_cat[1:5])

plot(association.rules_cat[1:10], method = "graph")

########### End R ############

