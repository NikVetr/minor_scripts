library(copula)
library(MASS)
library(psych)
library(rgl)
library(VineCopula)

set.seed(100)

m <- 3
n <- 2000
sigma <- matrix(c(1, 0.4, 0.2,
                  0.4, 1, -0.8,
                  0.2, -0.8, 1), 
                nrow=3)
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma)
cor(z,method='spearman')
pairs.panels(z)
u <- pnorm(z)
pairs.panels(u)
cor(u)
cor(z)
plot3d(u[,1],u[,2],u[,3],pch=20,col='navyblue')

x1 <- qgamma(u[,1],shape=2,scale=1)
x2 <- qbeta(u[,2],2,2)
x3 <- qt(u[,3],df=5)
plot3d(x1,x2,x3,pch=20,col='blue')

df <- cbind(x1,x2,x3)
pairs.panels(df)
cor(df,meth='spearman')

myCop <- normalCopula(param=c(0.4,0.2,-0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta", "t"),
              paramMargins=list(list(shape=2, scale=1),
                                list(shape1=2, shape2=2), 
                                list(df=5)) )

Z2 <- rMvdc(mvdc = myMvd, 2000)
colnames(Z2) <- c("x1", "x2", "x3")
pairs.panels(Z2)

cree <- read.csv('https://datascienceplus.com/wp-content/uploads/2015/10/cree_r.csv',header=F)$V2
yahoo <- read.csv('https://datascienceplus.com/wp-content/uploads/2015/10/yahoo_r.csv',header=F)$V2

plot(cree,yahoo,pch='.')
abline(lm(yahoo~cree),col='red',lwd=1)
cor(cree,yahoo,method='spearman')

u <- pobs(as.matrix(cbind(cree,yahoo)))[,1]
v <- pobs(as.matrix(cbind(cree,yahoo)))[,2]
pairs(cbind(u,v))
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula

t.cop <- tCopula(dim=2)
m <- pobs(as.matrix(cbind(cree,yahoo)))
fit <- fitCopula(t.cop,m,method='ml')
coef(fit)

rho <- coef(fit)[1]
df <- coef(fit)[2]
persp(tCopula(dim=2,rho,df=df),dCopula)

u <- rCopula(3965,tCopula(dim=2,rho,df=df))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')

cree_mu <- mean(cree)
cree_sd <- sd(cree)
yahoo_mu <- mean(yahoo)
yahoo_sd <- sd(yahoo)

hist(cree,breaks=80,main='Cree returns',freq=F,density=30,col='cyan',ylim=c(0,20),xlim=c(-0.2,0.3))
lines(seq(-0.5,0.5,0.01),dnorm(seq(-0.5,0.5,0.01),cree_mu,cree_sd),col='red',lwd=2)
legend('topright',c('Fitted normal'),col=c('red'),lwd=2)
hist(yahoo,breaks=80,main='Yahoo returns',density=30,col='cyan',freq=F,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(seq(-0.5,0.5,0.01),dnorm(seq(-0.5,0.5,0.01),yahoo_mu,yahoo_sd),col='red',lwd=2)
legend('topright',c('Fitted normal'),col=c('red'),lwd=2)

copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=df), margins=c("norm","norm"),
                    paramMargins=list(list(mean=cree_mu, sd=cree_sd),
                                      list(mean=yahoo_mu, sd=yahoo_sd)))
sim <- rMvdc(mvdc = copula_dist, 3965)

plot(cree,yahoo,main='Returns')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)


set.seed(4258)
copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=1), margins=c("norm","norm"),
                    paramMargins=list(list(mean=cree_mu, sd=cree_sd),
                                      list(mean=yahoo_mu, sd=yahoo_sd)))
sim <- rMvdc(mvdc = copula_dist, 3965)
plot(cree,yahoo,main='Returns')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated df=1'),col=c('black','red'),pch=21)
cor(sim[,1],sim[,2],method='spearman')

copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=8), margins=c("norm","norm"),
                    paramMargins=list(list(mean=cree_mu, sd=cree_sd),
                                      list(mean=yahoo_mu, sd=yahoo_sd)))
sim <- rMvdc(mvdc = copula_dist, 3965)
plot(cree,yahoo,main='Returns')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated df=8'),col=c('black','red'),pch=21)
cor(sim[,1],sim[,2],method='spearman')


#test for gene expression outliers?
myCop <- tCopula(param=c(rho.1 = -0.9), dim = 2, dispstr = "un")
# myCop <- tevCopula(param=c(rho.1 = -0.9), df = 0.1)
persp(myCop, dCopula)

myMvd <- mvdc(copula=myCop, margins=c("beta", "exp"),
              paramMargins=list(list(shape1=5, shape2=20),
                                list(rate=0.1)))
sim <- rMvdc(mvdc = myMvd, n = 10000)
sim[,1] <- sim[,1] / 2
plot(sim)
cor(sim)
