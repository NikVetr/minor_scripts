#functions
text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}

remove_bottom <- function(x, replacement){
  nobot <- gsub("g|j|p|q|y|,|\\(|\\)|Q", replacement, x)
  nobot <- gsub("\\_s*\\{[^\\)]+\\}", replacement, nobot) #underscore in brackets
  nobot <- gsub("_[a-z|0-9|A-Z]{1}", replacement, nobot) #underscore w/ one letter following
  nobot
}

remove_top <- function(x, replacement){
  notop <- gsub("\\^s*\\{[^\\)]+\\}", replacement, x)
  notop <- gsub("\\^[a-z|0-9|A-Z]{1}", replacement, notop)
  notop
}

remove_tb <- function(x, replacement){
  remove_top(remove_bottom(x, replacement), replacement)
}

text3 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, replacement = "a", ...){
  
  #convert text label to expression
  word_expression <- latex2exp::TeX(labels)
  
  #find general params
  strw <- strwidth(word_expression, cex = cex)
  strh <- strheight(word_expression, cex = cex)
  base_strh <- strheight(latex2exp::TeX("G"), cex = cex)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strh / 2, 0)
  
  #adjust in case of ledding
  nobot <- remove_bottom(labels, replacement)
  ebot <- strheight(latex2exp::TeX(nobot), cex = cex) - strh
  
  notop <- remove_top(labels, replacement)
  etop <- strheight(latex2exp::TeX(notop), cex = cex) - strh
  
  nobottop <- remove_tb(labels, replacement)
  ebottop <- strheight(latex2exp::TeX(nobottop), cex = cex) - strh
  
  #ugh this was obnoxious to figure out
  ebt_delta <- ebottop - (ebot + etop)
  adj_ledding <- ifelse(abs(ebt_delta) > 1E-6, 
                        ebot / 2 - (ebottop - ebot) / 2, 
                        ebot / 2 - etop / 2)
  adj_ledding <- adj_ledding - ifelse(base_strh > strh, (base_strh - strh) / 2, 0)
  
  #print the text itself
  text(x = adj_x, y = adj_y + adj_ledding, labels = word_expression, pos = NULL, cex = cex, col = col, ...)

  #draw a box around it if desired
  if(drect){
    rect(xleft = adj_x - strw / 2, 
         xright = adj_x + strw / 2, 
         ybottom = adj_y - strh / 2 + adj_ledding, 
         ytop = adj_y + strh / 2 + adj_ledding, border = col)
    abline(h=y - strheight(latex2exp::TeX("GIs"), cex = cex) / 2, lwd = 0.5)
  }
}

replacement <- "a"
plot(NULL, xlim = c(0,5), ylim = c(0,8), main = paste0("replacement = ", replacement))

cex = 10
xy1 <- c(0,0.125*cex)
xy2 <- c(0,0.5*cex)
xy3 <- c(0,0.825*cex)
# text(x = xy[1] + strwidth(lab, cex = cex) / 2, y = xy[2] + strheight(lab, cex = cex) / 2, labels = lab, cex = cex)
# rect(xleft = xy[1], xright = xy[1] + strwidth(lab, cex = cex), ybottom = xy[2], ytop = xy[2]+strheight(lab, cex = cex))

labs <- c("Testing", "$Testing_1$", "Testing$_1a^2$", "Testing$^2$", "Testing$_1^2$", ".")
# labs <- c("a", "$a^2$", "$a_1$", "$a_1^2$", "$a^2a_1$")
# labs <- c("a", "$a^.$", "$a_1$", "$a_1^.$", "$a^.a_1$")
# labs <- c("$_2$", "$^2$", "$_2^2")
cols <- adjustcolor(c("red", "green", "blue", "black", "orange", "lightgrey"), 0.5)

for(i in 1:length(labs)){
  text3(x = xy1[1], y = xy1[2], labels = labs[i], cex = cex, drect = T, pos = c(4), col = cols[i], replacement = replacement)
  text3(x = 3, y = 0.5+0.075*cex*i, labels = labs[i], cex = cex/5, drect = F, pos = c(4), col = cols[i], replacement = replacement)
  # text2(x = xy2[1], y = xy2[2], labels = latex2exp::TeX(labs[i]), cex = cex, drect = T, pos = c(4), col = cols[i])
  # text(x = xy3[1], y = xy3[2], labels = latex2exp::TeX(labs[i]), cex = cex, pos = c(4), col = cols[i])
}

points(xy1[1], xy1[2])



all_letters <- F
if(all_letters){
  cex = 1.5
  xy4s <- list(c(0,8), c(0,7), c(0,6))
  labs <- strsplit("\"!#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ^_`abcdefghijklmnopqrstuvwxyz|~", "")[[1]]
  labs <- c(labs, 
            paste0("$", c(LETTERS, letters), "^{2}$"), 
            paste0("$", c(LETTERS, letters), "_{2}$"), 
            paste0("$", c(LETTERS, letters), "^{2}_{1}$"))
  plotwidth <- diff(par("usr")[1:2]) / 1.25
  xy4adjx <- cumsum(c(0, strwidth(sapply(labs, latex2exp::TeX), cex = cex))) %% plotwidth
  xy4adjy <- floor(cumsum(c(0, strwidth(sapply(labs, latex2exp::TeX), cex = cex))) / plotwidth) / 2
  
  cols <- rep(1, length(labs))
  
  for(i in 1:length(labs)){
    text3(x = xy4[1] + xy4adjx[i], y = xy4[2] - xy4adjy[i], labels = labs[i], cex = cex, drect = T, pos = c(4), col = cols[i], replacement = replacement)
  }
  
}
