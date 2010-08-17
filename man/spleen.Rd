\name{spleen}
\alias{spleen}
\docType{data}
\title{Spleen Data set}
\description{
An example data set for paired data analysis.
A dataframe with the normalized Spleen data set as used in the original article.  
}
\usage{data(spleen)}
\format{
spleen is a dataframe. The first column contains gene names. The 6 columns A.1 to A.6 contain normalized mean log-intensity
and the 6 columns M.1 to M.6 correspond to normalized log-ratio
}
\references{P. Delmar, Robin, S., Tronik-Le Roux S. and Daudin J.-J. (2005) Mixture model on the variance for the differential analysis of gene expression data, JRSS series C, 54(1), 31:50}

\seealso{ \code{\link{DiffAnalysis}} }
\examples{
data(spleen)
# resSpleen <- DiffAnalysis(spleen,n=1,ind.array=1:6,name.A="A.",name.M="M.",
# varmixt.meth=TRUE,dyeswap=TRUE,indDS=c(2,4,6))
# Histogramm of the raw-pvalues
# hist(resSpleen$PValueVM)
}
\keyword{datasets}
