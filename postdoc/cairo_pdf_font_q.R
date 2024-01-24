cairo_pdf("test.pdf", width = 5, height = 5, family="Arial Unicode MS", pointsize = 18.5)
par(mar = c(0,0,0,0))
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim = c(0,1), ylim = c(0,1))
text(x = 0.5, y = 0.75, labels = "test 1", cex = 2, font = 3)
#oh wait it does work :/
dev.off()
