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

# Bibliography

- Croissant, Y., & Millo, G. (2018). Panel data econometrics with R. John Wiley & Sons.

- Kleiber, C., & Zeileis, A. (2008). Applied econometrics with R. Springer Science & Business Media.


# Acknowledgments

The author greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   
