#############################################
# Quantile Regression IV Panel Data Functions
# Harding and Lamarche, Economics Letters 104 (2009), 133-135 

##########################################################################################
# Function 1
##########################################################################################

rq.fit.fe <- function(X,y,w=1,taus=tau)
{
  
require(SparseM)
require(quantreg)
K <- length(w)
if(K != length(taus))
stop("length of w and taus must match")
X <- as.matrix(X)
p <- ncol(X)
n <- length(levels(as.factor(s)))
N <- length(y)
if( N != nrow(X))
 stop("dimensions of y,X must match")
D <- cbind(as(w,"matrix.diag.csr") %x% X)
y <- c(w %x% y)
a <- c((w*(1-taus)) %x% (t(X)%*%rep(1,N)))
rq.fit.sfn(D,y,rhs=a)

}

##########################################################################################
# Function 2
##########################################################################################

rq.fit.ivpanel <- function(d,exo,iv,y,s,tau)
{
  
# IV Quantile Regression with Fixed Effects
# s is an strata indicator
# d is the endogenous variable
# exo are the exogenous variables, no intercept.
# iv is the intrument
Z <- model.matrix(s~as.factor(s)-1)[,-1]
exo <- cbind(exo,Z)
X <- cbind(exo,1)
x <- cbind(d,X)
w <- cbind(iv,X)
ww <- t(w) %*% w
ww.inv <- ginv(as.matrix(ww))
wd <- t(w)%*%d
dhat <- w%*%ww.inv%*%wd
PSI <- cbind(dhat,X)
PSI1 <- cbind(d,X)
coef <- rq.fit.fe(PSI,y,tau=tau)$coef
resid <- y - PSI1%*%coef
mu1 <- mean(resid)
sigma1 <- var(resid)
c <- ((1-tau)*tau)/(dnorm(qnorm(tau,mu1,sqrt(sigma1)),mu1,sqrt(sigma1))^2)
PSIinv <- diag(length(coef))
PSIinv <- backsolve(qr(x)$qr[1:length(coef), 1:length(coef)], PSIinv)
PSIinv <- PSIinv %*% t(PSIinv)
vc1 <- c*PSIinv
std <- sqrt((length(y)/100)*diag(vc1))
alpha <- seq(coef[1]-2*std[1],coef[1]+2*std[1],by=std[1]/20)
z <- cbind(dhat,X)
betas <- matrix(NA,dim(z)[2],length(alpha))
g <- matrix(NA,length(alpha),1)
for (i in 1:length(alpha)){
ya <- y - alpha[i]*d
betas[,i] <- rq.fit.fe(z,ya,tau=tau)$coef
g[i,] <- max(svd(betas[1:dim(dhat)[2],i])$d)}
I <- which.min(g[,1])
param1 <- alpha[I]
est1 <- betas[(dim(dhat)[2]+1):dim(z)[2],I]
c(param1,est1)

}

##########################################################################################
# Function 3
##########################################################################################

rq.se.ivpanel <- function(bhat,d,exo,iv,y,s,tau)
{
  
n=length(y)
  
Z  <- model.matrix(s~as.factor(s)-1)[,-1]
X  <- cbind(exo,Z,1)
S0 <- cbind(iv,X)
D  <- cbind(d,X)
k  <- dim(D)[2]
vc <- matrix(0,k,k)
S  <- (1/n)*t(S0)%*%S0
  
res<- y - D%*%bhat
h  <- c(1.364 * ((2*sqrt(pi))^(-1/5))*sqrt(var(resid))*(n^(-1/5)))
J  <- (1/(n*h))*t(c(dnorm(resid/h)) %o% c(rep(1,dim(D)[2])) * D)%*%S0
vc <- (1/n)*(tau-tau^2)*ginv(as.matrix(J))%*%S%*%ginv(as.matrix(J))
  
rbind(bhat,(sqrt(diag(vc))))

}

##########################################################################################
##########################################################################################
# Libraries
#############################################

library(quantreg)
library(MASS)

#############################################
# A simple example:
#############################################

n   <- 10     # number of cross sectional units
m   <- 25     # number of time series units
rdu <- 0.5    # correlation parameter
tau <- 0.5    # quantile of interest
s <- rep(1:n,each=m)                         # strata 
Z <- model.matrix(s~as.factor(s)-1)[,-1]     # incidence matrix 

set.seed(25)
mu  <- rep(rchisq(n,3), each=m)
eps <- rchisq(length(s),3)
XZ  <- cbind(mu + epsilon,Z)
iv  <- rep(rnorm(n), each=m) + rnorm(length(s))
ief <- rep(rnorm(n), each=m)
v   <- rep(rnorm(n), each=m) + rnorm(length(s))
u   <- v*rdu + rnorm(length(s))*sqrt(1-rdu^2)
d   <- iv + v
y   <- ief + u

bhat  <- rq.fit.ivpanel(d,XZ[,1],iv,y,s,tau=tau)  
round(bhat,3)

bsehat <- rq.se.ivpanel(bhat,d,XZ[,1],iv,y,s,tau=tau)  
round(bsehat,3) -> bsehat

colnames(bsehat) <- c("endog", "exog", rep("ie",n-1), "int")
rownames(bsehat) <- c("beta", "se")
bsehat
##########################################################################################
