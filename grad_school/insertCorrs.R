ntraits <- 6
ind <- vector()
for (i in 1:(ntraits)){
  ind <- c(ind, i*(ntraits+1) - i*(i+1)/2)
}
ind

test <- rnorm(10, 1, 3)
indices <- seq(2, 10, by = 2)
insert <- rep(2, 5)
for(i in 1:length(insert)){
  test <- append(test, insert[i], after = indices[i]-1)
}
test

#apply to real data

strVec <- as.character(corrsUT)
ntraits <- ((8 * length(strVec) + 1) ^ 0.5 - 1) / 2 + 1

ind <- vector()
for (i in 1:(ntraits)){
  ind <- c(ind, i*(ntraits+1) - i*(i+1)/2)
}
ind

insert <- vector()
for (i in 1:length(ind)){
  insert[i] <- paste0("Corr[", i,"]")
}
insert

for(i in 1:length(insert)){
  strVec <- append(strVec, insert[i], after = ind[i]-1)
}
strVec
vec <- sapply(1:(length(strVec)-1), function (x) paste0(strVec[x], ","))
vec <- c(vec, strVec[length(strVec)])
vec[1] <- paste0("v(", vec[1])
vec[length(vec)] <- paste0(vec[length(vec)], ")")
cat(vec)

#find nearby correlation matrix
cor1 <- rbind(corrs, rep(0,length(corrs[,1])))
cor2 <- cbind(cor1, rep(0, length(cor1[,1])))
cor2[9,9] <- 1
cor2
