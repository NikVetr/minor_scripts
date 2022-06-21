library(dplyr)
library(microbenchmark)

x <- runif(1E8, 0, 1)
break_interval <- 0.001
breaks <- seq(0,1,by=break_interval)

fraction_smaller <- function(x, min, max, fineness) {
  a <- seq(min, max, by = fineness) 
  lx <- length(x) 
  la <- length(a) 
  lap <- 0 
  for (i in 1:la) { 
    b <- a[i] 
    below <- x < b 
    la[i] <- (sum(below) / lx) + lap 
    x <- x[!below] 
    lap <- la[i] 
  } 
  return(la) 
}

microbenchmark::microbenchmark(times = 10,
  (1-fraction_smaller(x = x, min = 0, max = 1, fineness = break_interval))*length(x),
  length(x) - c(0,cumsum(rle(sort(floor(x*1/break_interval)))$lengths)),
  sapply(breaks, function(i) sum(x > i)),
  ecdf(x),
  data.frame(cats = cut(x, breaks=c(breaks, Inf)), ordered_result=TRUE) %>% count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats)
)
