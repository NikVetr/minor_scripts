library(phytools)
library(phangorn)
library(phyloch)
install.packages("phyloch")


X <- read.dna(choose.files(), format = "sequential")
X <- as.phyDat(X)
Y <- read.dna(choose.files(), format = "sequential")
Y <- as.phyDat(Y)

objX <- pml(tree = rtree(n = length(X), tip.label = names(X), rooted = FALSE), 
           data = X)
objY <- pml(tree = rtree(n = length(Y), tip.label = names(Y), rooted = FALSE), 
           data = Y)
#objX <- update(objX, k=4)
#objY <- update(objY, k=4)

fit.X <- optim.pml(objX, optNNI = T, optEdge = T, optBf = T, optQ = T, optGamma = T, model = "GTR")
fit.Y <- optim.pml(objY, optNNI = T, optEdge = T, optBf = T, optQ = T, optGamma = T, model = "GTR")

par(mfrow = c(1,2))
plot(midpoint(fit.X$tree), main = "With Autapomorphies")
plot(midpoint(fit.Y$tree), main = "Without Autapomorphies")







X <- read.dna(choose.files(), format = "sequential")
X <- as.phyDat(X)
Y <- read.dna(choose.files(), format = "sequential")
Y <- as.phyDat(Y)

objX <- pml(tree = rtree(n = length(X), tip.label = names(X), rooted = FALSE), 
            data = X)
objY <- pml(tree = rtree(n = length(Y), tip.label = names(Y), rooted = FALSE), 
            data = Y)
#objX <- update(objX, k=4)
#objY <- update(objY, k=4)

fit.X <- optim.pml(objX, optNNI = T, optEdge = T, optBf = T, optQ = T, optGamma = T, model = "GTR")
fit.Y <- optim.pml(objY, optNNI = T, optEdge = T, optBf = T, optQ = T, optGamma = T, model = "GTR")

par(mfrow = c(1,2))
plot(midpoint(fit.X$tree), main = "With Autapomorphies")
plot(midpoint(fit.Y$tree), main = "Without Autapomorphies")













#Sankoff Parsimony

start <- rtree(n = 5, tip.label = names(B), br = NULL)
plotTree(start)
parsimony(start, B, method = "fitch")
mp.tree <- optim.parsimony(start, B, method = "fitch", all = TRUE)
plot(mp.tree, edge.width = 2)
mp.tree <- midpoint(mp.tree)  # you may have to do this twice
plotTree(mp.tree)

cost <- matrix(1, length(attr(primate.data, "levels")), length(attr(primate.data, 
                                                                    "levels")))
rownames(cost) <- colnames(cost) <- attr(primate.data, "levels")
diag(cost) <- 0
parsimony(mp.tree, primate.data, method = "sankoff", cost = cost)

transc <- 5
cost["a", "c"] <- cost["c", "a"] <- transc
cost["a", "t"] <- cost["t", "a"] <- transc
cost["c", "g"] <- cost["g", "c"] <- transc
cost["g", "t"] <- cost["t", "g"] <- transc
cost

mp.cost <- optim.parsimony(start, primate.data, method = "sankoff", cost = cost)

tree <- read.tree(text = "((A,B),(D,E));")
tree1 <- read.tree(text = "((A,B),((D,E),C));")
tree2 <- read.tree(text = "(((A,B),C),(D,E));")
par(mfrow = c(1,2))
plot(tree1)
plot(tree2)

tree1 <- pbtree(n=3)
tree2 <- rNNI(tree1)

A <- read.phyDat("C:\\Users\\Nikolai\\Dropbox\\autapomorphy_alignment_noC.dna", format = "phylip")
B <- read.phyDat("C:\\Users\\Nikolai\\Dropbox\\autapomorphy_alignment.dna", format = "phylip")

parsimony(tree1, B, method = "fitch")
parsimony(tree2, B, method = "fitch")

parsimony(tree1, B, method = "sankoff", cost = cost)
parsimony(tree2, B, method = "sankoff", cost = cost)

mp.cost <- optim.parsimony(tree2, A, method = "sankoff", cost = cost)
