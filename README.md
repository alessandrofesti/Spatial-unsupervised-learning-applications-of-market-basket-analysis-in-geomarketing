Spatial Association Rules For Geomarketing Purposes
================
Alessandro Festi
04 giugno 2019

<center>
Alessandro Festi
</center>
<br/>

<font size = '2'> Market basket analysis is a widely employed technique in marketing which provides suggestions of products to buy to a customer, given the past purchases made by the customers, basing upon the statistical methodology of association rules. This technique has been initially conceived to analyse transactions of products where the dataset is composed of sets of items purchased in different periods of time by N individuals, for instance formed from a large number of receipts collected in a point of sale of a large-scale retail trade. This chapter extends the application of market basket analysis to the geolocation points of N people in a specific space in order to discover associations among places that individuals have visited. When dealing with products, the goal is to link products while here it is linking locations. Similarly, one needs the sets of geolocation points for each individual considered in many different occasions. <br/> This code follows my article published by Springer </font>

``` r
library(geosphere)
library(rgeos)
library(shiny)
library(leaflet)
library(plyr) 
library(dplyr)
library(jsonlite)
library(arules)
library(MASS)
library(sp)
library(taRifx)
library(arulesViz)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(readxl)
```

An overview of the iperbole dataset, downloadable here: <http://dati.comune.bologna.it/node/640>

<center>
<img src="C:\Users\Alessandro\Desktop\OverviewIperboleDataset.jpg" width="90%" />
</center>
<br/> The R package 'reticulate' helps us integrating the Python language in the R environment. In the Iperbole dataset there are not latitude and longitude coordinates for the commercial activities in Bologna. We can get them using the geographic information we already have and geocoding it through the Mapbox API's using Python. Then the general analysis is implemented using R. <br/>

``` r
library(reticulate)
```

    ## Warning: package 'reticulate' was built under R version 3.4.4

``` r
use_python('C:\\Users\\afesti\\AppData\\Roaming\\Microsoft\\Windows\\Start Menu\\Programs\\Python 3.6')
```

``` python
import pandas as pd
import numpy as np 
from mapbox import Geocoder
import json

dataset = pd.read_csv("C:\\Users\\Alessandro\\elenco_esercizicommerciosedefissa_anno_2017.csv", sep = ';', header='infer', encoding='latin-1')

dataset['quartiere_settore'] = dataset.ESERCIZIO_VIA+'  '+dataset.ESERCIZIO_CIVICO+' '+dataset.QUARTIERE+' Bologna'

dataset['lat'] = float
dataset['lon'] = float

token = token
geocoder = Geocoder(access_token=token)

for i in range(len(dataset)):
    try:
        response = geocoder.forward(dataset.quartiere_settore[i])
        ale = response.content
        d = json.loads(ale)
        coordinates = d["features"][0]['geometry']['coordinates']
        dataset.iat[i,29] = coordinates[1]
        dataset.iat[i,30] = coordinates[0]
    except:
        dataset.iat[i,29] = np.nan
        dataset.iat[i,30] = np.nan

dataset.to_excel('geocoded.xlsx')
```

``` r
geocoded <- read_excel("C:\\Users\\Alessandro\\Desktop\\geocoded.xlsx")
geocoded <- dplyr::distinct(geocoded, lat, lon, .keep_all = TRUE)
```

<br/> Having obtained thir lat/lon coordinates, the position of the commercial activities in Bologna are plotted using Tableau <br/>
<center>
<img src="C:\Users\Alessandro\Desktop\MapBologna.png" width="90%" />
</center>
<br/>

Then the individual paths are simulated. In this case the number of simulated individuals is set to 50, each of them on 15 different occasions, geocoded 20 different times for each occasion. <br/>

``` r
n_individuals <- 50
n_paths <- 15
n_positions <- 20

points <- data.matrix(cbind(rnorm(n_individuals)/250 + 11.342220, rnorm(n_individuals)/250 + 44.493674)) 
correlation <- 0.7

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
```

    ##         lon      lat path ind
    ## 1  11.34135 44.48965    1   1
    ## 2  11.33560 44.48391    1   1
    ## 3  11.34226 44.49056    1   1
    ## 4  11.33359 44.48190    1   1
    ## 5  11.34292 44.49122    1   1
    ## 6  11.34101 44.48931    1   1
    ## 7  11.34370 44.49201    1   1
    ## 8  11.33761 44.48591    1   1
    ## 9  11.34011 44.48842    1   1
    ## 10 11.33975 44.48805    1   1

<br/> Then, in order to infer if a person was at a certain activity, the individual positions are matched with the positions of the geocoded commercial activities in Bologna on the basis of an arbitrary minimum distance(Euclidean distance).

``` r
path_v <- data.matrix(cbind(path_vis$lon,path_vis$lat))
places <- cbind(geocoded$lon, geocoded$lat)
places[is.na(places)] <- 0
min_dist <- 0.0005

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
```

    ##        lon      lat         dist path ind
    ## 1 11.34104 44.48999 0.0004532853    1   1
    ## 2       NA       NA 0.0025628323    1   1
    ## 3       NA       NA 0.0007576385    1   1
    ## 4       NA       NA 0.0034729058    1   1
    ## 5       NA       NA 0.0005741918    1   1
    ## 6       NA       NA 0.0005056982    1   1

<font size = '2'> Adjusting the dataset for the association rules discovery eliminating the non-matched observations </font>

``` r
result <- na.omit(result)
head(result)
```

    ##         lon      lat         dist path ind
    ## 1  11.34104 44.48999 0.0004532853    1   1
    ## 7  11.34412 44.49213 0.0004364515    1   1
    ## 12 11.34299 44.49179 0.0004977824    1   1
    ## 19 11.34095 44.48881 0.0004851239    1   1
    ## 20 11.34730 44.49525 0.0003898361    1   1
    ## 23 11.34184 44.48993 0.0002919664    2   1

Merging process to retrieve the index of the matched company

``` r
index <- dplyr::select(geocoded,lat,lon, IND)
index <- plyr::rename(index, c("lat" = "lat", "lon" = "lon", "IND" = "ID"))
final <- merge(x=index, y=result, by.x = c("lon","lat"), by.y = c("lon","lat"))
final$itemset_id <- paste(final$ind,final$path)
head(final,10)
```

    ##         lon      lat   ID         dist path ind itemset_id
    ## 1  11.32442 44.49021 3635 3.931014e-04    3  46       46 3
    ## 2  11.32442 44.49021 3635 2.088684e-04   13  46      46 13
    ## 3  11.32442 44.49021 3635 1.665054e-04   15  46      46 15
    ## 4  11.32452 44.49023 3663 1.088892e-04   14  46      46 14
    ## 5  11.32452 44.49023 3663 4.111399e-04    5  46       46 5
    ## 6  11.32452 44.49023 3663 7.369655e-05    4  46       46 4
    ## 7  11.32452 44.49023 3663 2.809274e-04    5  46       46 5
    ## 8  11.32527 44.49025 3654 4.796866e-04    4  46       46 4
    ## 9  11.32527 44.49025 3654 4.683413e-04   14  46      46 14
    ## 10 11.32527 44.49025 3654 4.827189e-04   15  46      46 15

Then one needs to transform the data into a transaction dataset

``` r
# Preparing for Association rules
transactionData <- ddply(final,c("itemset_id"),
                         function(final)paste(final$ID,
                                              collapse = ","))
```

The dataset is read again into R as a transaction object. The most frequent visited companies are displayed in the absolute frequency plot below

``` r
write.csv(transactionData,"C:\\Users\\Alessandro\\Desktop\\transaction_data.csv", quote = FALSE, row.names = TRUE)
tr <- read.transactions('C:\\Users\\Alessandro\\Desktop\\transaction_data.csv', format = 'basket', sep=',')
```

    ## Warning in asMethod(object): removing duplicated items in transactions

``` r
itemFrequencyPlot(tr,topN=20,type="relative",col=brewer.pal(8,'Pastel2'), main="Absolute Item Frequency Plot")
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png) <br/> Performing the Market Basket Analysis through the apriori algorithm

``` r
association.rules <- arules::apriori(tr, parameter = list(supp=0.01, conf=0.3,maxlen=10))
```

    ## Apriori
    ## 
    ## Parameter specification:
    ##  confidence minval smax arem  aval originalSupport maxtime support minlen
    ##         0.3    0.1    1 none FALSE            TRUE       5    0.01      1
    ##  maxlen target   ext
    ##      10  rules FALSE
    ## 
    ## Algorithmic control:
    ##  filter tree heap memopt load sort verbose
    ##     0.1 TRUE TRUE  FALSE TRUE    2    TRUE
    ## 
    ## Absolute minimum support count: 7 
    ## 
    ## set item appearances ...[0 item(s)] done [0.00s].
    ## set transactions ...[2046 item(s), 751 transaction(s)] done [0.00s].
    ## sorting and recoding items ... [281 item(s)] done [0.00s].
    ## creating transaction tree ... done [0.00s].
    ## checking subsets of size 1 2 3 4 done [0.00s].
    ## writing ... [1871 rule(s)] done [0.00s].
    ## creating S4 object  ... done [0.00s].

``` r
inspect(association.rules[1:10])
```

    ##      lhs       rhs    support    confidence lift     count
    ## [1]  {2144} => {2102} 0.01065246 0.8888889  33.37778  8   
    ## [2]  {2102} => {2144} 0.01065246 0.4000000  33.37778  8   
    ## [3]  {2144} => {2521} 0.01065246 0.8888889  17.56725  8   
    ## [4]  {1584} => {1856} 0.01065246 0.8000000  50.06667  8   
    ## [5]  {1856} => {1584} 0.01065246 0.6666667  50.06667  8   
    ## [6]  {5513} => {2521} 0.01331558 1.0000000  19.76316 10   
    ## [7]  {1696} => {2130} 0.01065246 0.8000000  26.12174  8   
    ## [8]  {2130} => {1696} 0.01065246 0.3478261  26.12174  8   
    ## [9]  {6855} => {2317} 0.01065246 0.6666667  17.26437  8   
    ## [10] {3259} => {3149} 0.01065246 0.6666667  18.54321  8

``` r
plot(association.rules, method = "two-key plot")
```

    ## To reduce overplotting, jitter is added! Use jitter = 0 to prevent jitter.

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
plot(association.rules, method = "graph")
```

    ## Warning: plot: Too many rules supplied. Only plotting the best 100 rules
    ## using 'support' (change control parameter max if needed)

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
plot(association.rules[1:20], method = "graph")
```

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png) <br/> The same approach is applied, rather than on company indexes, on their commercial sector

``` r
index_cat <- dplyr::select(geocoded,lat,lon, ATTIVITA_PREVALENTE_ESERCIZIO)
index_cat <- plyr::rename(index_cat, c("lat" = "lat", "lon" = "lon", "ATTIVITA_PREVALENTE_ESERCIZIO" = "Commercial type"))
final_cat <- merge(x=index_cat, y=result, by.x = c("lon","lat"), by.y = c("lon","lat"))
final_cat$itemset_id <- paste(final_cat$ind,final_cat$path)
final_cat <- na.omit(final_cat)
head(final_cat)
```

    ##         lon      lat                     Commercial type         dist path
    ## 11 11.32724 44.49036             prodotti ortofrutticoli 0.0003181488   15
    ## 12 11.32724 44.49036             prodotti ortofrutticoli 0.0002307451   10
    ## 58 11.32918 44.49045 Tabacco e altri generi di monopolio 0.0004924821    6
    ## 59 11.32918 44.49045 Tabacco e altri generi di monopolio 0.0004252687    7
    ## 60 11.32918 44.49045 Tabacco e altri generi di monopolio 0.0004423691    6
    ## 61 11.32918 44.49045 Tabacco e altri generi di monopolio 0.0004637482   15
    ##    ind itemset_id
    ## 11  47      47 15
    ## 12  47      47 10
    ## 58  36       36 6
    ## 59  36       36 7
    ## 60  36       36 6
    ## 61  48      48 15

``` r
transactionData_cat <- ddply(final_cat,c("itemset_id"),
                         function(final_cat)paste(final_cat$`Commercial type`,
                                              collapse = ","))
```

``` r
write.csv(transactionData_cat,"C:\\Users\\Alessandro\\Desktop\\transaction_data_cat.csv", quote = FALSE, row.names = TRUE)
tr_cat <- read.transactions('C:\\Users\\Alessandro\\Desktop\\transaction_data_cat.csv', format = 'basket', sep=',')
itemFrequencyPlot(tr_cat,topN=10,type="relative",col=brewer.pal(8,'Pastel2'), main="Absolute Item Frequency Plot")
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png) <br/> The positions of the commercial activities in Bologna are then plotted using Tableau according to their sector.

<center>
<img src="C:\Users\Alessandro\Desktop\MapBolognaByCategory.png" width="90%" />
</center>
<br/> Performing the Market Basket Analysis through the apriori algorithm on the sector

``` r
association.rules_cat <- arules::apriori(tr_cat, parameter = list(supp=0.03, conf=0.5,maxlen=5))
```

    ## Apriori
    ## 
    ## Parameter specification:
    ##  confidence minval smax arem  aval originalSupport maxtime support minlen
    ##         0.5    0.1    1 none FALSE            TRUE       5    0.03      1
    ##  maxlen target   ext
    ##       5  rules FALSE
    ## 
    ## Algorithmic control:
    ##  filter tree heap memopt load sort verbose
    ##     0.1 TRUE TRUE  FALSE TRUE    2    TRUE
    ## 
    ## Absolute minimum support count: 22 
    ## 
    ## set item appearances ...[0 item(s)] done [0.00s].
    ## set transactions ...[1658 item(s), 745 transaction(s)] done [0.00s].
    ## sorting and recoding items ... [34 item(s)] done [0.00s].
    ## creating transaction tree ... done [0.00s].
    ## checking subsets of size 1 2 3 done [0.00s].
    ## writing ... [15 rule(s)] done [0.00s].
    ## creating S4 object  ... done [0.00s].

``` r
inspect(association.rules_cat[1:5])
```

    ##     lhs                                rhs                            support confidence     lift count
    ## [1] {Biancheria per la casa}        => {Abbigliamento}             0.03489933  0.9285714 2.931295    26
    ## [2] {Mobili}                        => {Abbigliamento e accessori} 0.03892617  0.6041667 1.428902    29
    ## [3] {Cartoleria}                    => {Abbigliamento e accessori} 0.04697987  0.5645161 1.335125    35
    ## [4] {Pane - pasticceria - dolciumi} => {Abbigliamento}             0.04295302  0.5714286 1.803874    32
    ## [5] {Pane - pasticceria - dolciumi} => {Abbigliamento e accessori} 0.04429530  0.5892857 1.393707    33

``` r
plot(association.rules_cat, method = "graph")
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png) \`\`\`
