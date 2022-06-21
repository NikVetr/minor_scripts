cl <- c(3.692,3.681,3.7055,3.724,3.732,3.738,3.772,3.791,3.788)
mcl <- c(3.803,3.793,3.815,3.826,3.837,3.849,3.856,3.870,3.871)
scl <- c(3.900,3.905,3.914,3.915,3.924,3.931,3.936,3.940,3.946)
gpa <- cbind(cl = cl, mcl = mcl, scl = scl)
y <- 2013:2021
rownames(gpa) <- y

par(mar = c(6,6,4,2))
plot(1,1,xlim = c(2013, 2021.25), ylim = c(3.65,4), xaxt = "n", yaxt = "n", cex.lab = 2, cex.main = 2.5,
     xlab = "graduating year", ylab = "gpa threshold", main = "latin honors cutoffs @ vandy")
axis(1, at = y, labels = y, lwd = 2, cex.axis = 1.25)
axis(2, at = 72:80/20, labels = 72:80/20, lwd = 2, cex.axis = 1.25)
segments(x0 = y, x1 = y, y0 = 0, y1 = 100, lty = 3, col = "lightgrey", lwd = 0.25)
segments(x0 = 0, x1 = 1E4, y0 = 72:80/20, y1 = 72:80/20, lty = 3, col = "lightgrey", lwd = 0.25)
lines(y, cl, lwd = 3, col = "#DD8D29")
text(x = 2015, y = 3.7125, srt = 25, labels = "cum laude (top 25%)", font = 4, col = "#DD8D29")
lines(y, mcl, lwd = 3, col = "#46ACC8")
text(x = 2016.3, y = 3.8375, srt = 13.5, labels = "magna cum laude (top 13%)", font = 4, col = "#46ACC8")
lines(y, scl, lwd = 3, col = "#B40F20")
text(x = 2017.5, y = 3.9345, srt = 8, labels = "summa cum laude (top 5%)", font = 4, col = "#B40F20")
box(which = "plot", lwd = 4)

#add in school split
splits<-c(
"year","2019","2021","2020","
a&s","3.936","3.946","3.94","
a&s","3.856","3.871","3.87","
a&s","3.772","3.788","3.791","
bl","3.936","3.935","3.946","
bl","3.856","3.903","3.883","
bl","3.772","3.853","3.807","
eng","3.936","3.945","3.952","
eng","3.856","3.869","3.883","
eng","3.772","3.765","3.792","
pb","3.936","3.958","3.946","
pb","3.856","3.906","3.883","
pb","3.772","3.841","3.807")
splits <- gsub(pattern = "\n",replacement = "", splits)
splits <- matrix(splits, ncol = 4, byrow = T)
rownames(splits) <- splits[,1]
colnames(splits) <- splits[1,]
splits <- (splits[-1,-1])
class(splits) <- "numeric"
splits <- splits %*% matrix(c(1,0,0,0,0,1,0,1,0), nrow =3, byrow = F)
colnames(splits) <- 2019:2021
for(i in 1:nrow(splits)){
  lines(x = 2019:2021, y = splits[i,], lwd = 2, col = rep(c("#B40F20", "#46ACC8", "#DD8D29"), 4)[i])
}
text(rownames(splits), x = 2020.975, pos = 4, 
     col = rep(c("#B40F20", "#46ACC8", "#DD8D29"),3), cex = 1, y = splits[,"2021"])
