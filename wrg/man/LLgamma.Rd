\name{LLgamma}
\alias{LLgamma}

\title{Log Likelihood for gamma fitted data, failures and suspensions only}

\description{
 Determination of log-likelihood values for the basic distributions covered by R should be ubiquitous, but an implementation is provided here to enable likelihood ratio testing particularly for comparison of 3-parameter optimized models with their 2-parameter counterpart.
}

\usage{
LLgamma(x, s=NULL, alpha_g, beta_g)
}

\arguments{
  \item{x}{A vector of failure data. }
  \item{s}{An optional vector of suspension data.}
  \item{alpha_g}{The shape parameter from a gamma fit.}
  \item{beta_g}{The rate parameter from a gamma fit.}
}

\details{
 Function LLgamma can only handle datasets with failure and suspension data. For data including intervals \code{wblrLikelihood_g} can be used.
 This function is somewhat unique among likelihood functions in that it will scan the suspension data for
 removal of negative values. Such negative values are likely to be entered as data is provided
 for a 3-parameter fit by explicitly providing the original data with vector subtraction by the optimized t0.
 The primary intention for likelihood determination is to permit likelihood ratio testing for comparison
 of 3-parameter optimized models with their 2-parameter counterpart. 
}

\value{
Returns a log-likelihood value.
}

\references{
    Marie Laure Delignette-Muller, Christophe Dutang (2015). "fitdistrplus: An R Package for Fitting Distributions".
  Journal of Statistical Software, 64(4), 1-34. URL http://www.jstatsoft.org/v64/i04/.
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
LL<-LLgamma(failures, suspensions, 5.2 , 0.064)
}
