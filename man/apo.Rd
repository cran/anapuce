\name{apo}
\docType{data}
\alias{apo}
\title{Apo AI data }
\description{
An exemple data set for unpaired data analysis.
A dataframe with the normalized Apo data set as used in the original article.  
}
\usage{data(apo)}
\format{
apo is a dataframe. The first column contains gene names. The 8 columns cond1.1 to cond1.8 contain normalized normal mice measurements
and the 8 columns cond2.1 to cond2.8 correspond to normalized KO mice measurements
}
\references{M.J. Callow, S. Dudoit, E.L. Gong, T.P. Speed, and E.M. Rubin. Microarray expression profiling 
identifies genes with altered expression in hdl-deficien mice. Genome Res., 10(12) : 2022-9, 2000
 }
\seealso{ \code{\link{DiffAnalysis.unpaired}} }
\examples{
data(apo) 
# resApo <- DiffAnalysis.unpaired(apo,n=1,ind.array=1:8,varmixt.meth=TRUE)
# Histogramm of the raw-pvalues
# hist(resApo$PValueVM)
}
\keyword{datasets}