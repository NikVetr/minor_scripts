#PURPOSE is to get a sense of primate fossil abundance through time, efforts of discvovery, geological context, location, age

#discovery efforts
d <- read.csv("~/Downloads/PrimateFossilData.csv")
n <- length(d[,1])
uniq <- unique(d$accepted_name)
n_uniq <- length(uniq)
y_disc <- sapply(1:n_uniq, function(x) min(d$ref_pubyr[d$accepted_name == uniq[x]]))
discPerYear <- as.matrix(table(y_disc))
pubs <- rep(0, 2200)
pubs[as.numeric(rownames(discPerYear))] <- discPerYear[,1]
par(mfrow = c(1,2))
plot(1800:2019, pubs[1800:2019], type = "l", xlab = "Year of Publication", ylab = "Number of New Fossil Primate Species Published that Year", main = "Primate Fossil Discoveries")
plot(1800:2019, cumsum(pubs[1800:2019]), type = "l", xlab = "Year of Publication", ylab = "Number of Unique Fossil Primate Species Published in Total", main = "Primate Fossil Discoveries")

#geographical distribution
install.packages("rworldmap")
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, asp = 1)
points(d$lng, d$lat, col = 2, pch = 16, cex = 0.75)
uniq[startsWith(as.character(uniq), "Homo")]

#for epochs of the cenozoic, plot primate occurrence data
epoch_boundaries <- c(66, 56, 33.9, 23.03, 5.333, 2.58, 0.012) # in MA, from GSA 2018 https://www.geosociety.org/GSA/Education_Careers/Geologic_Time_Scale/GSA/timescale/home.aspx
epoch_names <- c("Paleocene", "Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene", "Holocene")
par(mar=c(0,0,0,0))
par(mfrow = c(4,2))
for(i in 1:length(epoch_names)){
  print(i)
  epoch_name <- epoch_names[i]
  plot(newmap, asp = 1)
  text(labels = epoch_name, x = 120, y = -80)
  
  #subset points and plot
  #if really fancy, use geographic envelopes. But we are not really fancy so points will do for now
  start <- epoch_boundaries[i]
  end <- epoch_boundaries[i+1]
  sub_d <- d[(d$max_ma <= start & d$max_ma >= end) | (d$min_ma <= start & d$min_ma >= end),] 
  points(sub_d$lng, sub_d$lat, col = 2, pch = 16, cex = 0.75)
}
