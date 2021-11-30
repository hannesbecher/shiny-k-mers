# Descriptions of the data sets.

#' K-mer spectrum of \emph{Euphrasia arctica} individual E028
#'
#' A k-mer spectrum if allotetraploid \emph{E. arctica}.
#'
#' @format A \code{spectrum} object containing a name and data frame:
#' \describe{
#'   \item{`name`}{The name 'E. arctica, E028'}
#'   \item{`k`}{The k-mer length of 27}
#'   \item{`data`}{A dataframe with coulmns \describe{
#'   \item{`mult`}{K-mer multiplicity}
#'   \item{`count`}{The number of different k-mer at a given multiplicity}
#'   }
#'   }
#' }
#' @source Becher et al. (2020) \url{https://doi.org/10.1016/j.xplc.2020.100105}
#'
#' @examples plot(E028, log="xy")
#' plot(E028, xlim=c(0, 200), ylim=c(0, 10000000))
#' tetmer(E028)
"E028"



#' K-mer spectrum of \emph{Euphrasia anglica} individual E030
#'
#' A k-mer spectrum of diploid \emph{E. anglica}.
#'
#' @format A \code{spectrum} object containing a name and data frame:
#' \describe{
#'   \item{`name`}{The name 'E. anglica, E030'}
#'   \item{`k`}{The k-mer length of 21}
#'   \item{`data`}{A dataframe with coulmns \describe{
#'   \item{`mult`}{K-mer multiplicity}
#'   \item{`count`}{The number of different k-mer at a given multiplicity}
#'   }
#'   }
#' }
#' @source Source!
#'
#' @examples plot(E030, log="xy")
#' plot(E030, xlim=c(0, 200), ylim=c(0, 10000000))
#' tetmer(E030)
"E030"

