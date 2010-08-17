\name{LBI}
\alias{LBI}
\title{ Label Bias Index }
\description{
This function calculates the LBI (Label Bias Index) for two self-self hybridization microarrays.
}
\usage{
LBI(infile, name.M = "M.norm", ind.array = 1:2, graph = TRUE,
graphout = "FigM1M2")
}
\arguments{
  \item{infile      }{ input data file with normalized log-ratio base 2.}
  \item{name.M      }{ character string containing a \link{regular expression} contained in the columnames corresponding to the M values.}
  \item{ind.array   }{ the indices of arrays to analyze. }
  \item{graph       }{ logical indicating wether to perform graphics. }
  \item{graphout    }{ name of the graphic output. }
}

\references{ Martin-Magniette M.L, Aubert J, Cabannes E, Daudin J.J (2005) Evaluation of the gene-specific dye bias in cDNA microarray experiments, \emph{Bioinformatics}, \bold{21(9)}, 1995-2000.  }
\author{ J. Aubert }

\keyword{ programming }
