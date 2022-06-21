gdp <- read.csv("data/gdp_pc.csv")
geo <- read.csv("data/Country_List_ISO_3166_Codes_Latitude_Longitude.csv")
conts <- read.csv("data/countriesContinents.csv")

gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

deg2rad <- function(deg) return(deg*pi/180)

gdp <- gdp[!is.na(gdp$gdp_percap),]
compConts <- intersect(geo$Alpha.3.code, gdp$iso3c)
d <- matrix(0, nrow = length(compConts), ncol = 2)
colnames(d) <- c("distance_from_ethiopia_km", "ln_gdp_percap_ppp_int$")
countries <- rep("country", length(compConts))
continents <- rep("continent", length(compConts))
waypoints = T

for(i in 1:nrow(d)){
  c1 <- "ETH"
  c2 <- compConts[i]
  
  c1latlong <- c(deg2rad(geo$Latitude..average.[geo$Alpha.3.code == c1]), deg2rad(geo$Longitude..average.[geo$Alpha.3.code == c1]))
  c2latlong <- c(deg2rad(geo$Latitude..average.[geo$Alpha.3.code == c2]), deg2rad(geo$Longitude..average.[geo$Alpha.3.code == c2]))
  c2cont <- conts$region[conts$alpha.3 == c2]
  d[i,2] <- log(gdp$gdp_percap[gdp$iso3c == c2])
  countries[i] <- as.character(gdp$country[gdp$iso3c == c2])
  continents[i] <- as.character(c2cont)
  
  if(waypoints){
    AnadyrRussia <- c(64, 177)
    CairoEgypt <- c(30, 31)
    IstanbulTurkey <- c(41, 28)
    PhnomPenhCambodia <- c(11, 104) 
    PrinceRupertCanada <- c(54, -130)
    if(c2cont == "Africa"){
      d[i,1] <- gcd.hf(c1latlong[2], c1latlong[1], c2latlong[2], c2latlong[1])
    }
    if(c2cont == "Europe"){
      d1 <- gcd.hf(c1latlong[2], c1latlong[1], CairoEgypt[2], CairoEgypt[1])
      d2 <- gcd.hf(CairoEgypt[2], CairoEgypt[1], IstanbulTurkey[2], IstanbulTurkey[1])
      d3 <- gcd.hf(IstanbulTurkey[2], IstanbulTurkey[1], c2latlong[2], c2latlong[1])
      d[i,1] <- sum(d1,d2,d3)
    }
    if(c2cont == "Asia"){
      d1 <- gcd.hf(c1latlong[2], c1latlong[1], CairoEgypt[2], CairoEgypt[1])
      d2 <- gcd.hf(CairoEgypt[2], CairoEgypt[1], c2latlong[2], c2latlong[1])
      d[i,1] <- sum(d1,d2)    
    }
    if(c2cont == "Oceania"){
      d1 <- gcd.hf(c1latlong[2], c1latlong[1], CairoEgypt[2], CairoEgypt[1])
      d2 <- gcd.hf(CairoEgypt[2], CairoEgypt[1], PhnomPenhCambodia[2], PhnomPenhCambodia[1])
      d3 <- gcd.hf(PhnomPenhCambodia[2], PhnomPenhCambodia[1], c2latlong[2], c2latlong[1])
      d[i,1] <- sum(d1,d2,d3)
    }
    if(c2cont == "Americas"){
      d1 <- gcd.hf(c1latlong[2], c1latlong[1], CairoEgypt[2], CairoEgypt[1])
      d2 <- gcd.hf(CairoEgypt[2], CairoEgypt[1], AnadyrRussia[2], AnadyrRussia[1])
      d3 <- gcd.hf(AnadyrRussia[2], AnadyrRussia[1], PrinceRupertCanada[2], PrinceRupertCanada[1])
      d4 <- gcd.hf(PrinceRupertCanada[2], PrinceRupertCanada[1], c2latlong[2], c2latlong[1])
      d[i,1] <- sum(d1,d2,d3,d4)
    }
  } else {
    d[i,1] <- gcd.hf(c1latlong[2], c1latlong[1], c2latlong[2], c2latlong[1])
  }
}

library(RColorBrewer)
continentColors <- unique(continents)
continentColors <- cbind(continentColors, brewer.pal(length(continentColors), "Dark2"))
plot(d, col = sapply(1:length(continents), function(x) continentColors[which(continents[x] == continentColors[,1]),2]))
legend(35000, 12, legend = continentColors[,1], fill = continentColors[,2])
cor(d)

