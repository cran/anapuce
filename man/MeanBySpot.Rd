\name{MeanBySpot}
\alias{MeanBySpot}
\title{ Mean By Spot }
\description{
  Computes the mean on the different values of spots grouping by \emph{by.var}.
}
\usage{
MeanBySpot(fileIN, n = 3, name.A = "A", name.M = "M.norm", 
by.var = "ID", na.rm=TRUE)
}
\arguments{
  \item{fileIN  }{ input dataframe. }
  \item{n       }{ number of identificant columns. }
  \item{name.A  }{ character string containing a \link{regular expression} contained in the columnames corresponding to the A values. }
  \item{name.M  }{ character string containing a \link{regular expression} contained in the columnames corresponding to the M values. }
  \item{by.var  }{ name of the grouping variable. }
  \item{na.rm   }{ a logical value indicating whether NA values should be stripped before the computation proceeds}
}
\value{
  A R dataframe
}
\author{ J. Aubert }
\seealso{ \code{\link{DiffAnalysis}}, \code{\link{DiffAnalysis.unpaired}},\code{\link{mean}} }
\keyword{ manip }
