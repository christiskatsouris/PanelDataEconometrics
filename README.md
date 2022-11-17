# PanelDataEconometrics

In this teaching page, we focus on estimation and inference methodologies for panel data models using R. The data structure can involve either time-series or non time-series related data. Specifically for panel time series regressions we initially operate under the assumption of stationarity. In particular, an existing helpful R package in the [plm](https://cran.r-project.org/web/packages/plm/vignettes/A_plmPackage.html) package.


# Estimation Examples

### Example 1: Panel Data Regression Model 

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
