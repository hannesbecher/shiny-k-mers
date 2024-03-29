% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E030}
\alias{E030}
\title{K-mer spectrum of \emph{Euphrasia anglica} individual E030}
\format{
A \code{spectrum} object containing a name and data frame:
\describe{
\item{\code{name}}{The name 'E. anglica, E030'}
\item{\code{k}}{The k-mer length of 21}
\item{\code{data}}{A dataframe with coulmns \describe{
\item{\code{mult}}{K-mer multiplicity}
\item{\code{count}}{The number of different k-mer at a given multiplicity}
}
}
}
}
\source{
Source!
}
\usage{
E030
}
\description{
A k-mer spectrum of diploid \emph{E. anglica}.
}
\examples{
plot(E030, log="xy")
plot(E030, xlim=c(0, 200), ylim=c(0, 10000000))
tetmer(E030)
}
\keyword{datasets}
