\name{MEV}

\alias{Regression in model with errors in the variables}

\title{Regression in model with errors}

\usage{
MEV( x, y, lambda_x)
}

\description{
MEV is used to adjust models with errors in the variables.
It can be used to perform regressions on variables that have measurement errors.
}


\examples{
MEV(seq(1,100,length.out = 100), cumsum(runif(100,1,100)), 2)
}
