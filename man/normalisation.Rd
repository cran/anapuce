\name{normalisation}
\alias{normalisation}
\title{ Function to normalize microarray data }
\description{
  Performs the data normalization using a general loess and 
  a block effect correction via the subtraction of the median by block if there is a Block column.
}
\usage{
normalisation(fileIN = NULL, Red = "F635.Median", Green = "F532.Median",
n = 3, flag = -100, graph = TRUE, filter.function = filterByDefault, 
filter.fic = NULL, filter.var = NULL, sep.write = "\t", dec.write = ".",
header = TRUE, sep = "\t", skip = 0, ...)
}
\arguments{
  \item{fileIN          }{ on optional regular expression contained in names of the files to analyze}
  \item{Red             }{ name of the column containing the Cy5 values}
  \item{Green           }{ name of the column containing the Cy3 values }
  \item{n               }{ number of identificant columns }
  \item{flag            }{ list of flags values to delete }
  \item{graph           }{ logical indicating wether to perform graphics }
  \item{filter.function }{ function to perform before beginning the normalization }
  \item{filter.fic      }{ argument passed to filter.funtion - file containing the values of the variable named \emph{filter.var} to delete }
  \item{filter.var      }{ argument passed to filter.funtion }
  \item{sep.write       }{ the field separator string to use in the \link{write.table} function }
  \item{dec.write       }{ the string to use for decimal points in the \link{write.table} function }
  \item{header          }{ a logical value indicating whether the file contains the names of the variables as its first line. Used in the \link{read.table} function }
  \item{sep             }{ the field separator string to use in the \link{read.table} function }
  \item{skip            }{ integer: the number of lines of the data file to skip before beginning to read data. }
  \item{\dots           }{ Further arguments to be passed to read.table. }
}
\details{
  By default \emph{fileIN=NULL} and the files to analyze have to be choosen among all the files in the working directory. 
  If fileIN is not NULL, only file names which match the \link{regular expression} fileIN will be proposed to the normalization step.
}
\value{
  A file \emph{resNorm.txt} containing the normalized data and if \emph{graph=TRUE} a file .pdf per array with graphics
}
\references{ Yang, Y., Dudoit, S., Luu, P., Lin, D., Peng, V., Ngai, J. and Speed, T., (2002), Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation, \emph{Nucleic Acids Research} \bold{30}.
 }
\author{ J. Aubert }

\seealso{ \code{\link{filterByDefault}},\code{\link{read.table}}}
\keyword{ programming }

