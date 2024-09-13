cols <- colors(distinct = T)
crgb <- t(col2rgb(cols))
lens <- apply(crgb, 1, function(x) sqrt(sum(x^2)))
crgbs <- cbind(crgb * 1/lens, c(scale(lens)))
dists <- dist(crgbs)
hc <- hclust(dists, method = "complete") 
n <- 50
clusters <- cutree(hc, k = n)
crgb_split <- split(data.frame(crgb), clusters)
separated_cols <- do.call(rbind, lapply(crgb_split, function(x){
  average_dist <- sqrt(colSums(as.matrix(dist(x))^2)) / (nrow(x)-1)
  most_connected_col <- x[which.min(average_dist),]
  return(most_connected_col)
}))
selected_cols <- apply(separated_cols[order(cmdscale(dist(separated_cols), k = 1)),], 1, function(x) 
  rgb(x[1], x[2], x[3], maxColorValue = 255))
plot(1:n, col = selected_cols, pch = 19, cex = 2)
