\name{filterByDefault}
\alias{filterByDefault}
\title{ Filter applied on data before normalisation }
\description{
  Apply a filter on input data by deleting spots with values for filter.var are mentionned in \emph{filter.fic}
                               by deleting spots with flags equal to \emph{flag0} (-100 by default).
}
\usage{
filterByDefault(fileIN, flag0, m, filter.fic = filter.fic, filter.var = filter.var, sep = sep, sep.write = sep.write, dec.write = dec.write, ...)
}

\arguments{
  \item{fileIN      }{ input datafame. }
  \item{flag0       }{ value(s) of flags to delete. }
  \item{m           }{ number of identificant columns. }
  \item{filter.fic  }{ name of the file containing the values of the variable named filter.var to delete. }
  \item{filter.var  }{ name of the variable concerned with deletion. }
  \item{sep         }{ the field separator string to use in the \link{read.table} function. }
  \item{sep.write   }{ the field separator string to use in the \link{write.table} function. }
  \item{dec.write   }{ the string to use for decimal points in the \link{write.table} function. }
  \item{\dots       }{ Further arguments to be passed to read.table. }
}
\value{
  A subset of the input dataframe
}
\author{J. Aubert}

\seealso{ \code{\link{normalisation}}}

\keyword{ manip }

