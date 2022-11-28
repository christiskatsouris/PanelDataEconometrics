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

### Example 3: Quantile Panel Data Regressions (MC Simulation)

Reference to the following example, is the paper of Galvao, A., Lamarche, C., and Lima, L. R. (JASA, 2013).

```R

# Initialization of parameters
#######################################################
beta1 <- 10             # First slope
beta2 <- -2             # Second slope
beta  <- c(beta1,beta2) # Parameters of interest
gamma <- 0.5          # Model with heteroskedastic errors
Ci    <- -0.95        # Censoring point
tau   <- 0.5          # Quantile of interest
const <- 0.05         # Initial value for d in 3-step
deltan <- const       # Initial value for c of 3-step
Rep    <- 1000        # Number of repetitions
#######################################################

datagen <- function(n,t,beta1,beta2,gamma,Ci)
{# begin-of-function
 	
s <- rep(1:n,rep(t,n))
X <-array(0,c(n*t,2))
	
for (i in 1:(n*t))
{# begin-for-loop
   signal=0
   while (signal==0)
   {# begin-while-loop
     x <- rnorm(2)
     if (max(abs(x))<=2) 
       {X[i,]<-x;signal=1}
   }# end-while-loop
}# end-for-loop
 	
X1  <- X[,1]
X2  <- X[,2]
X1s <- X1^2
X2s <- X2^2
X   <- cbind(X1,X2,X1s,X2s)

eta.aux <- rnorm(n)
medX    <- (X[,1]+X[,2])/2
eta     <- array(0,c(n,1))
eta[1]  <- eta.aux[1] + ( 1/sqrt(1) )*sum( medX[1:t] )
  
for (k in 1:(n-1))
{
  eta[k+1]<-eta.aux[k+1]+(1/sqrt(1))*sum(medX[(k*t+1):((k+1)*t)])
}
 
eta   <- rep(eta.aux,rep(t,n))
u     <- rnorm(n*t)
ystar <- eta + beta1*X1+ beta2*X2 + (1+gamma*(X1+X2+X1s+X2s))*u
y     <- replace(ystar, ystar < Ci, Ci)
delta <- 1-((y==Ci)*1)
  
return(cbind(y,X,s,delta,ystar))
}# end-of-function 
``` 

## Key References: 

$\textbf{[1]}$ Galvao, A. F., Juhl, T., Montes-Rojas, G., & Olmo, J. (2018). Testing slope homogeneity in quantile regression panel data with an application to the cross-section of stock returns. Journal of Financial Econometrics, 16(2), 211-243.

$\textbf{[2]}$ Härdle, W. K., Wang, W., & Yu, L. (2016). Tenet: Tail-event driven network risk. Journal of Econometrics, 192(2), 499-513.

$\textbf{[3]}$ Galvao, A. F., Lamarche, C., & Lima, L. R. (2013). Estimation of censored quantile regression for panel data with fixed effects. Journal of the American Statistical Association, 108(503), 1075-1089.

$\textbf{[4]}$ Kato, K., Galvao Jr, A. F., & Montes-Rojas, G. V. (2012). Asymptotics for panel quantile regression models with individual effects. Journal of Econometrics, 170(1), 76-91.

$\textbf{[5]}$ Lamarche, C. (2010). Robust penalized quantile regression estimation for panel data. Journal of Econometrics, 157(2), 396-408.

$\textbf{[6]}$ Harding, M., & Lamarche, C. (2009). A quantile regression approach for estimating panel data models using instrumental variables. Economics Letters, 104(3), 133-135.


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
