\name{wblrLoglike}
\alias{wblrLoglike}

\title{Log likelihood for Weibull, Lognormal and Gamma fitted data including intervals}

\description{This function generates the Log likelihood for Wiebull, lognormal and gamma distributions from life(-time)
 data containing fixed time occurrences, and suspensions (progressive right-censored data), discoveries (left-censored data),
 and interval censored data. This function will overwrite existing wblrLoglike in WeibullR.
}
\usage{
wblrLoglike(par, x, dist="weibull", sign=1, tz=0 )
}

\arguments{
  \item{par}{A vector of parameters taken in the same order as the underlying R distribution.
  That is, for weibull the order of parameters is (shape, scale).  For lognormal (meanlog, sdlog).
  (Unsuitable parameters cannot be pre-tested, but may cause the function to fail with NaN, or inf values.)
  }
  \item{x}{A dataframe such as generated by mleframe with column names 'left', 'right' and optionally 'qty'.
    Exact failure data (occurrences) have same time entered in both 'left' and 'right' columns.
	Suspension data has last known life[time] entered in 'left' column and -1 entered in 'right' column.
	The left(early) interval bound for left-censored data must be entered as zero. (NA is not accepted).
  }
  \item{dist}{A string defining a distribution to be fit. Implemented distributions are "weibull" (default), "lognormal" and "gamma".
  (Only 2-parameter models are accepted.)
  }
  \item{sign}{An integer of 1 or -1 assigning a multiplier for the returned value. A value of -1 results in a minimization
  function suitable for use as an argument to optim.
  }
  \item{tz}{A scalar vector with the translation parameter, t0, to be applied to data.
  }
}

\value{
  Returns a scalar vector holding the value of log-likelihood.  Unsuitable parameters may cause generation of NaN.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 

  Marie Laure Delignette-Muller, Christophe Dutang (2015). "fitdistrplus: An R Package for Fitting Distributions".
  Journal of Statistical Software, 64(4), 1-34. URL http://www.jstatsoft.org/v64/i04/.
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
weibull_loglike<-wblrLoglike(c( 2.26, 4900), mleframe(failures,suspensions))
}
