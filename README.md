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


# Acknowledgments

The author greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   
