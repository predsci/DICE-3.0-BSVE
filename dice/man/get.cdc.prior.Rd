% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_fxns.R
\name{get.cdc.prior}
\alias{get.cdc.prior}
\title{Retrieve Prior Parameters for a CDC Run}
\usage{
get.cdc.prior(mydata = NULL)
}
\arguments{
\item{mydata}{A dataframe with all the available for this \pkg{DICE} run}
}
\value{
mydata dataframe with additional information about the prior
}
\description{
\code{get.cdc.prior} is a wrapper that calls the \code{\link{get.prior}} function
}
\examples{
mydata <- get.cdc.prior(mydata = mydata)

}
