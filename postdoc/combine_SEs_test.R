var2 <- function(x, bessel = T) sum((x - mean(x))^2) / (length(x) - ifelse(bessel, 1, 0))

foo <- function(n1, n2){
  x1 <- rnorm(n1)
  x2 <- rnorm(n2)
  dx <- (mean(x1) - mean(x2))
  tx <- dx / (sqrt(var2(x1,T) / n1 + var2(x2,T) / n2))
  pval <- (1 - pt(abs(tx), n1 + n2 - 2)) * 2
  
  return(c(my = pval, tt = t.test(x1, x2)$p.value))
}

n1 <- 10
n2 <- 15

out <- t(replicate(100, foo(n1, n2)))
plot(out); abline(0,1)
max(abs(out[,1] - out[,2]))
