Functions for time series exploration, modelling and forecasting for R: tsutils package
=======
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tsutils?color=blue)](https://CRAN.R-project.org/package=tsutils)
[![Downloads](http://cranlogs.r-pkg.org/badges/tsutils?color=blue)](https://CRAN.R-project.org/package=tsutils)


Development repository for the tsutils package for R.

Stable version can be found at: https://cran.r-project.org/package=tsutils

## Installing

To install the development version use:
```{r}
if (!require("devtools")){install.packages("devtools")}
devtools::install_github("trnnick/tsutils")
```
Otherwise, install the stable version from CRAN:
```{r}
install.packages("tsutils")
```
## Functionality

The **tsutils** package provides functions to support various aspects of time series and forecasting modelling. In particular this package includes: (i) tests and visualisations that can help the modeller explore time series components and perform decomposition; (ii) modelling shortcuts, such as functions to construct lagmatrices and seasonal dummy variables of various forms; (iii) an implementation of the Theta method; (iv) tools to facilitate the design of the forecasting process, such as ABC-XYZ analyses; and (v) "quality of life" tools, such as treating time series for trailing and leading values.

#### Time series exploration:
* `cmav`: centred moving average.
* `coxstuart`: Cox-Stuart test for location/dispersion.
* `decomp`: classical time series decomposition.
* `seasplot`: construct seasonal plots.
* `trendtest`: test a time series for trend.

#### Time series modelling:
* `getOptK`: optimal temporal aggregation level for AR(1), MA(1), ARMA(1,1).
* `lagmatrix`: create leads/lags of variable.
* `residout`: construct control chart of residuals.
* `seasdummy`: create seasonal dummies.
* `theta`: Theta method.

#### Forecasting process modelling:
* `abc`: ABC analysis.
* `xyz`: XYZ analysis.
* `abcxyz`: ABC-XYZ analyses visualisation.

#### Quality of life:
* `geomean`: geometric mean.
* `lambdaseq`: generate sequence of lambda for LASSO regression.
* `leadtrail`: remove leading/training zeros/NAs.
* `wins`: winsorisation, including vectorised versions `colWins` and `rowWins`.

#### Time series data:
* `referrals`: A&E monthly referrals.

## Authors & contributors

* Nikolaos Kourentzes - (http://nikolaos.kourentzes.com/)
* Ivan Svetunkov - (https://forecasting.svetunkov.ru/)
* Oliver Schaer - (https://www.lancaster.ac.uk/lums/people/oliver-schaer)

## References
References are provided where necessary at the help file of each function. The overall modelling philosophy is reflected in:

Ord K., Fildes R., Kourentzes N. (2017) [Principles of Business Forecasting, 2e](http://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/). Wessex Press Publishing Co.

## License

This project is licensed under the GPL3 License

_Happy forecasting!_
