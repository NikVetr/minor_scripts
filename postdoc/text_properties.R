
plot(NULL, xlim = c(0,5), ylim = c(0,8))
cex = 10
xy <- c(2,3)
# text(x = xy[1] + strwidth(lab, cex = cex) / 2, y = xy[2] + strheight(lab, cex = cex) / 2, labels = lab, cex = cex)
# rect(xleft = xy[1], xright = xy[1] + strwidth(lab, cex = cex), ybottom = xy[2], ytop = xy[2]+strheight(lab, cex = cex))

text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, ...){
  
  #convert text label to expression
  word_expression <- latex2exp::TeX(labels)
  
  #find general params
  strw <- strwidth(word_expression, cex = cex)
  strh <- strheight(word_expression, cex = cex)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strh / 2, 0)
  
  #adjust in case of ledding
  ebot <- strheight(latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", labels)), cex = cex) - strh
  etop <- strheight(latex2exp::TeX(gsub("\\^", "a", labels)), cex = cex) - strh
  ebottop <- strheight(latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u|\\^", "a", labels)), cex = cex) - strh
  
  #this only works if prefix and suffix are not adjacent, otherwise things get moved
  adj_ledding <- ifelse(abs(ebottop - (ebot + etop)) > 1E-6, 
                        (ebot + etop - ebottop) * 2 / 10, 
                        ebot / 2 - etop / 2)
  
  #ugh this was obnoxious to figure out
  adj_ledding <- ifelse(abs(ebottop - (ebot + etop)) > 1E-6, 
                        ebot / 2 - (ebottop - ebot) / 2, 
                        ebot / 2 - etop / 2)
  
  #print the text itself
  text(x = adj_x, y = adj_y + adj_ledding, labels = word_expression, pos = NULL, cex = cex, col = col, ...)
  
  #draw a box around it if desired
  if(drect){
    rect(xleft = adj_x - strw / 2, 
         xright = adj_x + strw / 2, 
         ybottom = adj_y - strh / 2 + adj_ledding, 
         ytop = adj_y + strh / 2 + adj_ledding, border = col)
  }
}

lab <- "Sest"
text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(4), col = adjustcolor(1, 0.5))

lab <- "Sest$_1$"
text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(4), col = adjustcolor("green", 0.5))

lab <- "Sest$_1a^2$"
text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(4), col = adjustcolor("blue", 0.5))

lab <- "Sest$^2$"
text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(4), col = adjustcolor("red", 0.5))

lab <- "Sest$_1^2$"
text2(x = xy[1], y = xy[2], labels = lab, cex = cex, drect = T, pos = c(4), col = adjustcolor("orange", 0.5))
# 

points(xy[1], xy[2])
