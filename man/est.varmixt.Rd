\name{est.varmixt}
\alias{est.varmixt}
\title{Variance Mixture Estimation}
\description{
Performs variance mixture analysis
}
\usage{
est.varmixt(VAR, Kmax, dfreedom)
}
\arguments{
  \item{VAR     }{ vector of estimated variance. }
  \item{Kmax    }{ maximal number of variance components. }
  \item{dfreedom}{ degrees of freedom of the estimated variance. }
}

\value{
 a LIST with the following components :
  \item{BIC.crit}{value of the BIC criterion}
  \item{p.i     }{the probability of each variance component}
  \item{vars    }{variances}
  \item{loglike }{value of the criterion based on Loglikelihood}
  \item{nmixt   }{number of variance components}
  \item{tau     }{the matrix of posterior probability that a gene belongs to each variance component. One row per gene, one column per variance component. }
  \item{VM2     }{the variance attributed to each gene according to the MAP rule - Delmar et al. (2005) \emph{JRSS}}
  \item{VM      }{the variance attributed to each gene taking into account the tau values - Delmar et al. (2005) \emph{Bioinformatics}} 
}
\references{ 
Delmar P, Robin S, Le Roux D, Daudin J.J (2005), Mixture model on the variance for the differential analysis of gene expression, \emph{JRSS series C}, \bold{54:1}, 31-50.

Delmar P, Robin S, Daudin J.J (2005), VarMixt: efficient variance modelling for the differential analysis of replicated gene expression data, \emph{Bioinformatics}, \bold{21(4)}, 502-8.
}
\author{ M-L Martin-Magniette and J. Aubert }

\seealso{ \code{\link{DiffAnalysis}},\code{\link{DiffAnalysis.unpaired}} }
\keyword{ programming }
