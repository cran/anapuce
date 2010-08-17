\name{FileChoice}
\alias{FileChoice}
\title{ Function to choose the files to analyze }
\description{
  This function sets the working directory to the directory choosen by the user,
  and allows the user to choose the files to analyze inside this directory.
}
\usage{
FileChoice(pattern = NULL)
}
\arguments{
  \item{pattern     }{ an optional \link{regular expression}. Only file names which match the regular expression will be returned. }
}

\value{
  A vector containing the names of the files to analyze
}
\author{ J. Aubert }

\seealso{ \code{\link{normalisation}} }

\keyword{ internal }
