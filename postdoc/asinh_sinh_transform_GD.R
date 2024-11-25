x <- -1000:1000/100
asx <- asinh(x)
sx <- sinh(x)

ps <- -10:10/10
avelines <- lapply(ps, function(p) sx * p + asx * (1-p))
plot(x, asx, type = "l", lwd = 2, col = 3)
lines(x, sx, lwd = 2, col = 3)
for(i in 1:length(ps)){
  lines(x, avelines[[i]])
}
abline(0,1,col=2,lty=2)
