% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_fxns.R
\name{psi.lnlike}
\alias{psi.lnlike}
\title{Calculate the Likelihood of a Model}
\usage{
psi.lnlike(mod, dat)
}
\arguments{
\item{mod}{Numeric time series with the predicted number of cases}

\item{dat}{Integer time series with the observed number of cases}
}
\value{
myLLK  The total likelihood of all the data points
}
\description{
Calculate the total likelihood given two time-series: the observed one and predicted (dat and mod respectively).
}
\examples{
psi.lnlike(mod=model_rtn, dat = mydata$model$epi)
psi.lnlike(mod=model_rtn, dat = mydata$fit$epi[,iregion])

}
