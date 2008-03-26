\name{DiffAnalysis.unpaired}
\alias{DiffAnalysis.unpaired}
\title{ Differential Analysis for unpaired data }
\description{
  Performs the differential analysis from normalised unpaired data according to different ways of variance modelling and
  computes lists of differentially expressed genes according different multiple testing procedures.
}
\usage{
DiffAnalysis.unpaired(fileIN = "resNorm.txt", n = 3, cond1 = "cond1.", cond2 = "cond2.", fileOUT = "ListOfGenes.txt", fileDelete = "GenesOutOfAnalysis.txt", procs = c("bonferroni", "BH"), alpha = c(0.05, 0.05), fileID = NULL, function.trt = NULL, by.var = "ID", varmixt.meth = FALSE, header = TRUE, sep = "\t", sep.write = "\t", dec.write = ".", ...)
}

\arguments{
  \item{fileIN      }{ normalized data set. }
  \item{n           }{ number of identificant columns. }
  \item{cond1       }{ a \link{regular expression} corresponding to the first condition. }
  \item{cond2       }{ a \link{regular expression} corresponding to the second condition. }
  \item{fileOUT     }{ output data file. }
  \item{fileDelete  }{ output data file containing the list of withdrawn genes. }
  \item{procs       }{ adjusting p-values procedures to be used. }
  \item{alpha       }{ values of the first type error to be used in the different procedures. }
  \item{fileID      }{ file giving information about genes. }
  \item{function.trt}{ function to be applied before differential analysis. }
  \item{by.var      }{ argument passed to function.trt. }
  \item{varmixt.meth}{ logical indicating whether to perform the variance modelling proposed by Delmar et al. (2005). }
  \item{header      }{ a logical value indicating whether the file contains the names of the variables as its first line. Used in the \link{read.table} function. }
  \item{sep         }{ the field separator string to use in the \link{read.table} function. }
  \item{sep.write   }{ the field separator string to use in the \link{write.table} function. }
  \item{dec.write   }{ the string to use for decimal points in the \link{write.table} function. }
  \item{\dots       }{ Further arguments to be passed to read.table. }
}
\details{
  This function performs a differential analysis in the gene-specific and homoscedastic cases.
  
  If \emph{varmixt.meth = TRUE}, the method used is the method proposed by Delmar et al. (2005), (the same as coded in the vm.analysis function from the varmixt package.)

}
\value{
  By default these following files are created
  
  - one list of differentially expressed genes by procedure

  - The list of withdrawn genes

  - The list of all the genes
 
  An R object is returned if function return value is assigned.

}
\references{ Delmar, P., Robin, S. and Daudin, J.J., (2005), VarMixt: efficient variance modelling for the differential analysis of replicated gene expression data, \emph{Bioinformatics}, \bold{21},(4), 502--8
   
   Dudoit, S., Yang, Y. H., Callow, M. J. and Speed, T.P., (2002), Statistical methods for identifying differentially expressed genes in replicated cdna microarray experiments, \emph{Statisticia Sinica}, \bold{12}, 111--139
 }
\author{ J. Aubert }

\seealso{ \code{\link{p.adjust}},  \code{\link{MeanBySpot}},  \code{\link{DiffAnalysis}} }
\examples{

}

\keyword{ programming }

