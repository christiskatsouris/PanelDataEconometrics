# Panel Data Econometrics

In this teaching page, we focus on estimation and inference methodologies for panel data models using R. The data structure can involve either time-series or non time-series related data. Specifically for panel time series regressions we initially operate under the assumption of stationarity. In particular, an existing helpful R package in the [plm](https://cran.r-project.org/web/packages/plm/vignettes/A_plmPackage.html) package.


# Estimation and Inference Examples

### Example 1: Panel Data Regression Model with plm R package

```R

# Definition of Design matrix 
X <- model.matrix(x)

# Condition 
if (nrow(X) <= ncol(X)) stop("insufficient number of observations")

y  <- pmodel.response(x)
r  <- lm(y ~ X - 1, model = FALSE)
nc <- colnames(model.frame(r)$X)
names(r$coefficients) <- nc

``` 

### Example 2: Univariate Cross-Section Quantile Regression Model (Tenet Dataset)

```R

dataset <- read.csv("100_firms_returns_and_macro_2015-04-15.csv", header = TRUE)

returns <- as.matrix(dataset[, 2:101])    
macro   <- as.matrix(dataset[, 102:108])

firm.data <- read.csv(file = "Bal_sheet.csv")
firm.data <- as.matrix(firm.data)

# Take a subset of only 15 firms
returns <- returns[ , 1:15]

for (k in 1:15) 
{
  firm.data.new <- firm.data[ , 1:(4 * k)]
}  

ncol <- NCOL(firm.data.new)
return.t   <- as.matrix( returns[2:314 ,1] )
return.lag <- as.matrix( returns[1:313, 1] )

## Step 1: Estimate Univariate Quantile Regression Models with covariates y_{t-1} and firm characteristics

k <- 1
firm.data.new <- firm.data[2:314 , (4 * k - 3):(4 * k)]
model.quantile         <- rq( return.t  ~ return.lag + firm.data.new, tau = 0.05 )
model.quantile.summary <- summary( model.quantile, se = "boot", bsmethod= "xy" )

> model.quantile.summary

Call: rq(formula = return.t ~ return.lag + firm.data.new, tau = 0.05)

tau: [1] 0.05

Coefficients:
                    Value    Std. Error t value  Pr(>|t|)
(Intercept)          0.25567  1.10568    0.23124  0.81728
return.lag           0.10764  0.17216    0.62524  0.53228
firm.data.newLEV.1  -0.05732  0.01896   -3.02379  0.00271
firm.data.newMM.1    0.77781  0.54157    1.43622  0.15196
firm.data.newSIZE.1  0.04860  0.02733    1.77802  0.07639
firm.data.newMTB.1   0.01098  0.08260    0.13289  0.89437

``` 

### Example 3: Quantile Panel Data Regressions 

```R

#############################################
# Quantile Regression IV Panel Data Functions
# Harding and Lamarche, Economics Letters 104 (2009), 133-135 
##########################################################################
# Functions
##########################################################################
rq.fit.fe <- function(X,y,w=1,taus=tau){
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
rq.fit.sfn(D,y,rhs=a)}
#############################################
rq.fit.ivpanel <- function(d,exo,iv,y,s,tau){
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
c(param1,est1)}
#############################################
rq.se.ivpanel <- function(bhat,d,exo,iv,y,s,tau){
n=length(y)
Z <- model.matrix(s~as.factor(s)-1)[,-1]
X <- cbind(exo,Z,1)
S0 <- cbind(iv,X)
D <- cbind(d,X)
k=dim(D)[2]
vc <- matrix(0,k,k)
S <- (1/n)*t(S0)%*%S0
resid <- y - D%*%bhat
h <- c(1.364 * ((2*sqrt(pi))^(-1/5))*sqrt(var(resid))*(n^(-1/5)))
J = (1/(n*h))*t(c(dnorm(resid/h)) %o% c(rep(1,dim(D)[2])) * D)%*%S0
vc = (1/n)*(tau-tau^2)*ginv(as.matrix(J))%*%S%*%ginv(as.matrix(J))
rbind(bhat,(sqrt(diag(vc))))}
##########################################################################
##########################################################################
# Libraries
library(quantreg)
library(MASS)
#############################################
# A simple example:
#
n <- 10                                      # number of cross sectional units
m <- 25                                      # number of time series units
rdu = 0.5                                    # correlation parameter
tau = 0.5                                    # quantile of interest
s <- rep(1:n,each=m)                         # strata 
Z <- model.matrix(s~as.factor(s)-1)[,-1]     # incidence matrix 
set.seed(25)
mu <- rep(rchisq(n,3), each=m)
epsilon <- rchisq(length(s),3)
XZ <- cbind(mu + epsilon,Z)
iv <- rep(rnorm(n), each=m) + rnorm(length(s))
ief <- rep(rnorm(n), each=m)
v <- rep(rnorm(n), each=m) + rnorm(length(s))
u <- v*rdu + rnorm(length(s))*sqrt(1-rdu^2)
d <- iv + v
y <- ief + u
rq.fit.ivpanel(d,XZ[,1],iv,y,s,tau=tau) -> bhat
round(bhat,3)
rq.se.ivpanel(bhat,d,XZ[,1],iv,y,s,tau=tau) -> bsehat
round(bsehat,3) -> bsehat
colnames(bsehat) <- c("endog", "exog", rep("ie",n-1), "int")
rownames(bsehat) <- c("beta", "se")
bsehat

``` 

## Key References: 

$\textbf{[1]}$ Kato, K., Galvao Jr, A. F., & Montes-Rojas, G. V. (2012). Asymptotics for panel quantile regression models with individual effects. Journal of Econometrics, 170(1), 76-91.

$\textbf{[2]}$ Lamarche, C. (2010). Robust penalized quantile regression estimation for panel data. Journal of Econometrics, 157(2), 396-408.

$\textbf{[3]}$ Harding, M., & Lamarche, C. (2009). A quantile regression approach for estimating panel data models using instrumental variables. Economics Letters, 104(3), 133-135.


# Bibliography

$\textbf{[1]}$ Croissant, Y., & Millo, G. (2018). Panel data econometrics with R. John Wiley & Sons.

$\textbf{[2]}$ Pesaran, M. H. (2015). Time series and panel data econometrics. Oxford University Press.

$\textbf{[3]}$ Kleiber, C., & Zeileis, A. (2008). Applied econometrics with R. Springer Science & Business Media.

$\textbf{[4]}$ Arellano, M. (2003). Panel data econometrics. OUP Oxford.


# Acknowledgments

The author greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   
