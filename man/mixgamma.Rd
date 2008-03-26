\name{mixgamma}
\alias{mixgamma}
\title{ Gamma mixture }
\description{
This function provides an estimation of a mixture of gamma distributions.
The code is inspired from the varmixt package.
}
\usage{
mixgamma(VAR, dfreedom, var.init, pi.init, nmixt, stop.crit , display = TRUE, niter.max = 50000, criterion = criterion)
}
\arguments{
  \item{VAR         }{ vector of estimated variance. }
  \item{dfreedom    }{ degrees of freedom of the estimated variance. }
  \item{var.init    }{ vector of initial variances for the mixture. }
  \item{pi.init     }{ vector of initial prior probability for the mixture. }
  \item{nmixt       }{ integer: number of components in the mixture model. }
  \item{stop.crit   }{ the stopping relative precision limit for stopping EM algorithm. }
  \item{display     }{ boolean: should the result of model fitting be displayed on the screen ? }
  \item{niter.max   }{ integer : maximum number of iterations of the EM algorithm. }
  \item{criterion   }{ criterion for the convergence of the EM algorithm : either "likelihood" 
                 for criterion base on loglikelihood either "parameter". }
}
\details{
  
}
\value{
  a LIST with the following components :
  \item{BIC.crit}{value of the BIC criterion}
  \item{p.i     }{the probability of each variance component}
  \item{vars    }{variances}
  \item{loglike }{value of the criterion base on Loglikelihood}
  \item{nmixt   }{number of variance components}
  \item{tau     }{the matrix of posterior probability that a gene belongs to each variance component. One row per gene, one column per variance component. }
  \item{VM2     }{the variance attributed to each gene according to the MAP rule - Delmar et al. (2005) \emph{JRSS}}
  \item{VM      }{the variance attributed to each gene taking into account the tau values - Delmar et al. (2005) \emph{Bioinformatics}} 
}
\references{ }
\author{ Paul Delmar and M-L Martin-Magniette}
\note{ 
}
\seealso{ \code{\link{est.varmixt}}}
\examples{
  }
\keyword{ internal }
