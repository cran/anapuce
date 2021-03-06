\name{LocalFDR}
\alias{LocalFDR}
\title{ Local FDR }
\description{Estimation of local false discovery rates
}
\usage{
LocalFDR(dataf = dataf, graph = TRUE, method = NULL, lambda0 = 0.5, 
smoothing = "1", thres = c(0.01, 0.05, 0.1, 0.2), mm = c(3, 5, 15, NA))
}

\arguments{
  \item{dataf       }{ input data file with two columns (1 = gene name, 2= ordered raw pvalues). }
  \item{graph       }{ if TRUE a pdf graphic is created.}
  \item{method      }{ method for estimating m0. This must be one of the strings "adaptive", "conservative", "bootstap", "smoother" or NULL. }
  \item{lambda0     }{ value used in calculating m0. }
  \item{smoothing   }{ "1" for the initial published method, "2" for the PAVA method
                      isotonic (monotonely increasing nonparametric) least squares regression - see P. Broberg (2005). }
  \item{thres       }{ threshold defining intervals used for the movering average smoothing. }
  \item{mm          }{ parameter defining intervals used for the movering average smoothing. }
}

\value{
 - if graph = TRUE, a file of graphics named LocalFDRGraph.pdf

 - A data file (LocalFDRFile.txt) with 4 columns : gene name, raw pvalues and two columns corresponding to smoothed FDR values.
}
\references{ Aubert J, Bar-Hen A, Daudin J.J, Robin S (2004) Determination of the differentially expressed genes in microarrays experiments using local FDR, \emph{BMC Bioinformatics}, \bold{5:125}.

Aubert J , Bar-Hen A, Daudin J.J, Robin S (2005) Correction: Determination of the differentially expressed genes in microarray experiments using local FDR, \emph{BMC Bioinformatics}, \bold{6:42}.

 Per Broberg (2005) A comparative review of estimates of the proportion unchanged genes and the false discovery rate, \emph{BMC Bioinformatics} \bold{6:199} }

\author{ J. Aubert }

\keyword{ programming }

