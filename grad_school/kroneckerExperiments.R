dim1 <- 15
dim2 <- 20
A <- matrix(rnorm(dim1^2), dim1)
B <- matrix(rnorm(dim1^2), dim1)
C <- matrix(rnorm(dim2^2), dim2)
D <- matrix(rnorm(dim2^2), dim2)
ABAt <- A%*%B%*%t(A)
CDCt <- C%*%D%*%t(C)
any(abs(kronecker(ABAt, CDCt) - kronecker(A, C) %*% kronecker(B, D) %*% t(kronecker(A, C))) > 0.000001)

det(kronecker(ADAt, BEBt))
det(ADAt)^dim(BEBt)[1] * det(BEBt)^dim(ADAt)[1]

A <- matrix(rnorm(dim1^2), dim1)
B <- matrix(rnorm(dim1^2), dim1)
X <- rnorm(dim1*dim1)
kronecker(A,B)%*%X

A%*%rep(1,15)
rbind(cbind(A,A),cbind(A,A))%*%rep(1,30)

A <- matrix(rnorm(dim1^2), dim1)
X <- rnorm(dim1)
Am <- cbind(A, A%*%X)
Xm <- c(X, -1)
Am%*%Xm
rbind(cbind(Am,Am),cbind(Am,Am))%*%rep(Xm,2)
rbind(cbind(3*Am,1.1*Am),cbind(17*Am,0*Am))%*%rep(Xm,2)

dim1 <- 10
A <- matrix(rnorm(dim1^2), dim1)
X <- rnorm(dim1)
A <- cbind(A, A%*%X)
X <- c(X, -1)
A%*%X #evaluates to 0, as expected
dim2 <- 5
B <- matrix(rnorm(dim2^2), dim2)
kronecker(B, A)%*%rep(X,dim2) #evaluates to 0, as also expected
