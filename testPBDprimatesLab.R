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
par(mfrow = c(1,1))
plot(1800:2019, pubs[1800:2019], type = "l", xlab = "Year of Publication", ylab = "Number of New Fossil Primate Species Published that Year", main = "Primate Fossil Discoveries")
plot(1800:2019, cumsum(pubs[1800:2019]), type = "l", xlab = "Year of Publication", ylab = "Number of Unique Fossil Primate Species Published in Total", main = "Primate Fossil Discoveries")

#geographical distribution
install.packages("rworldmap")
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, asp = 1)
points(d$lng, d$lat, col = 2, pch = 16, cex = 0.75)
uniq[startsWith(as.character(uniq), "Homo")]
