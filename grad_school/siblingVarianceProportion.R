n <- 1e5
varA <- 45
varB <- 14
varC <- 14
a = rnorm(n=n,m=0,sd=varA^0.5)
b = rnorm(n=n,m=0,sd=varB^0.5)
c = rnorm(n=n,m=0,sd=varC^0.5)
# create a series of sib pairs for which a common factor (a) ought to explain ~50% of variance in "phenotype"
sib1 = a+b
sib2 = a+c
sib1 = sib2 - c + b
cor(sib1, sib2)
cor(a, sib2)^2
lm(sib1 ~ sib2)$coefficients[2]

# d <- - c + b
# var(d) #should be varC + varB
# 
# # so 1 - varD / var(A) + var(B) = cor(sib1, sib2)^2
# 1 - (varC + varB)/(varA + varB)
# cor(sib1, sib2)^2
# #close but not quite