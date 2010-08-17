\name{rowVars}
\alias{rowVars}
\title{ Form row variances}
\description{
Form row variances for numeric arrays. 
}
\usage{
rowVars(x, ...)
}
\arguments{
  \item{x}{ an array of two or more dimensions, containing numeric, complex, integer or logical values, or a numeric data frame.
}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
This function is equivalent tu use of apply with FUN = var with MARGIN = 1, but is a lot faster.
}
\value{
A numeric or complex array of suitable size, or a vector if the result is one-dimensional.
}
\seealso{ \code{\link{rowSums}}, \code{\link{rowMeans}}, \code{\link{apply}}  }
\keyword{  internal}
