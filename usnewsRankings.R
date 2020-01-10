setwd("~/Documents/USNews")

#read in data + key for different naming conventions betw. d1 vs. (d2, d3)
d1 <- read.csv("usnews_83-07.csv", check.names=FALSE); rownames(d1) <- d1[,1]; d1 <- d1[-58,-1]; colnames(d1)
d2 <- read.csv("usnews_08-15.csv", check.names=FALSE); rownames(d2) <- d2[,1]; d2 <- d2[,-1]; d2 <- d2[,1:5]; colnames(d2)
d3 <- read.csv("usnews_13-20.csv", check.names=FALSE); rownames(d3) <- d3[,1]; d3 <- d3[,-1]; colnames(d3)
matcher <- read.csv("usnews_namematch.csv", header = F)

#add in schools not in key / d1
addInNewSchools <- T
if(addInNewSchools){
  notInD1 <- setdiff(unique(rownames(d2), rownames(d3)), matcher[,2])
  matcher <- as.data.frame(cbind(c(as.character(matcher$V1), notInD1), c(as.character(matcher$V2), notInD1)))
}

#make some names nicer
matcher <- as.matrix(matcher)
matcher[matcher[,1] == "Northeastern",1] <- "Northeastern University"
matcher[matcher[,1] == "Boston Univ",1] <- "Boston University"
matcher[matcher[,1] == "Georgia",1] <- "University of Georgia"
matcher <- as.data.frame(matcher)

#compile main data frame
d <- t(sapply(1:length(matcher[,1]), function(x) c(d1[rownames(d1) == matcher[x,1],], d2[rownames(d2) == matcher[x,2],], d3[rownames(d3) == matcher[x,2],])))
rownames(d) <- matcher[,1]
d <- as.data.frame(d)

#get positions for names in right margin (since US News allows ties)
name_positions <- orig_name_positions <- sort(unlist(d[,"2020"]))
for(i in 2:length(name_positions)){if(name_positions[i] <= name_positions[i-1]){name_positions[i] <- name_positions[i-1] + 1}}
name_positions[-(1:15)] <- name_positions[-(1:15)] - 1

#set plotting window properties
par(mar = c(6,6,0,15))
par(xpd=TRUE)

#specify colors for rainbow effect or not
cols <- rep(1, length(d[,1]))

#subset schools for easier viewing
subsetD <- T
maxRank <- 50
schoolSize <- 1.2
if(subsetD){
  d <- d[as.numeric(which(d$`2020` < (maxRank + 1))),]
}

#throw away data for small handful of low-ranked schools
d <- d[sapply(1:length(d[,1]), function(x) which(rownames(d) == labels(name_positions)[x])),]
d[d > 65] <- NA

#draw plots
dir.create("plots")
for(j in nrow(d):1){
  print(paste0(j, ": ", rownames(d)[j]))
  png(filename = paste0("plots/", j, "_", rownames(d)[j], ".png"), width = 1200, height = 1200)
  par(mar = c(3,5,0,20.5))
  par(xpd=TRUE)
  plot(1,1, bty="n", xlim = c(1984,2020), ylim = c(1,65), type = "n", xaxt = "n", yaxt="n", xlab = "", ylab = "")
  title("U.S. News National University Rankings (1983 - 2020)", line = -4, adj = 0.05, cex.main = 3)
  title(ylab="Rank", cex.lab=2)
  title(xlab = "Year", line=1.5, cex.lab=2)
  axis(side = 2, at = c(1:13*5-4), cex.axis=1.5)
  axis(side = 1, at = 1983:2020, pos = 0, cex.axis = 1.25)
    for(i in (1:length(d[,1]))[-j]){
      years <- as.numeric(colnames(d))
      rank <- d[i,]
      name <- rownames(d)[i]
      color <- cols[i]
      line_width <- 1
      lines(years[!is.na(rank)], rank[!is.na(rank)], col = color, lwd = line_width)
      textStart <- 2020
      if(sum(orig_name_positions[labels(orig_name_positions) == name] == orig_name_positions) > 1){
        textStart <- 2022
        segments(x0 = 2020, y0 = orig_name_positions[labels(orig_name_positions) == name], 
                 x1 = 2022, y1 = name_positions[labels(name_positions) == name], lty = 3, lwd = 1.5)
      }
      text(name, x = textStart, y = name_positions[labels(name_positions) == name], pos = 4, col = color, cex = schoolSize)
    }
  color = 2
  line_width <- 3
  years <- as.numeric(colnames(d))
  rank <- d[j,]
  name <- rownames(d)[j]
  lines(years[!is.na(rank)], rank[!is.na(rank)], col = color, lwd = line_width)
  textStart <- 2020
  if(sum(orig_name_positions[labels(orig_name_positions) == name] == orig_name_positions) > 1){
    textStart <- 2022
    segments(x0 = 2020, y0 = orig_name_positions[labels(orig_name_positions) == name], col = color,
             x1 = 2022, y1 = name_positions[labels(name_positions) == name], lty = 3, lwd = 1.5)
  }
  text(name, x = textStart, y = name_positions[labels(name_positions) == name], pos = 4, col = color, cex = schoolSize)
  dev.off()
}

#put together gif
library(magick)
library(gtools)
library(purrr)

rev(mixedsort(list.files(path = "plots/", pattern = "*.png", full.names = T))) %>% 
  map(image_read) %>% 
  image_join() %>%  
  image_animate(fps=1) %>% 
  image_write(path = "schools.gif", format = "gif") 

