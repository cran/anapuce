\name{Intensity.Norm}
\alias{Intensity.Norm}
\title{ Function to compute normalized intensity for each channel }
\description{
  Computes normalized intensity for each channel and each array from the file containing normalized data.
  }
\usage{
Intensity.Norm(fileIN = "resNorm.txt", n = 3, ind.array = NULL, name.A = "A", name.M = "M.norm", sep = "\t", center = FALSE, log.transf = TRUE, ...)
}
\arguments{
  \item{fileIN      }{ input dataframe. }
  \item{n           }{ number of identificant columns. }
  \item{ind.array   }{ indices of arrays to use. }
  \item{name.A      }{ character string containing a \link{regular expression} contained in the columnames corresponding to the A values. }
  \item{name.M      }{ character string containing a \link{regular expression} contained in the columnames corresponding to the M values. }
  \item{sep         }{ the field separator string to use in the \link{read.table} function.}
  \item{center      }{ logical indicating wether data will be centered by array. }
  \item{log.transf  }{ logical indicating wether data will stay log-transformed. }
  \item{\dots}{ Further arguments to be passed to read.table. }
}
\details{
}
\value{
  A R dataframe containing the n first columns of the input dataframe and the "normalized" intensity 
  for the red and green channels of each chosen array.
}
\references{Thorne, N.P. (2004) Single-channel normalisation and analysis of twocolour cDNA microarray data, PhD thesis }
\author{ J. Aubert }
\seealso{ \code{\link{normalisation}}}

\keyword{ programming }

