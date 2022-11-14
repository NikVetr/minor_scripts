
plot(NULL, xlim = c(0,5), ylim = c(0,8))
lab <- "Sest"
cex = 10
xy <- c(3,3)
# text(x = xy[1] + strwidth(lab, cex = cex) / 2, y = xy[2] + strheight(lab, cex = cex) / 2, labels = lab, cex = cex)
# rect(xleft = xy[1], xright = xy[1] + strwidth(lab, cex = cex), ybottom = xy[2], ytop = xy[2]+strheight(lab, cex = cex))

text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(lab, cex = cex) / 2, 
         xright = adj_x + strwidth(lab, cex = cex) / 2, 
         ybottom = adj_y - strheight(lab, cex = cex) / 2, 
         ytop = adj_y + strheight(lab, cex = cex) / 2)
  }
}

text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(2,4))
points(xy[1], xy[2])
