#' Plot the confidence intervals of Gini coefficient
#'
#' Computes and plot the normalized Gini coefficient (G*) (Raffinetti et al.
#' 2015) and its confidence intervals of the residual frequency distributions
#' of the  desired statistics: \code{"geoD"} (Geodesic Distances),
#' \code{"paco"} (PACo) or \code{"paraF"} (M03)
#' (either all or only one of them).
#'
#' @param LF_1 Vector of statistics produced with
#'        \code{\link[=max_cong]{max_cong()}} or
#'        \code{\link[=max_incong]{max_incong()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"}.
#'
#' @param M01 Matrix produced with
#'        \code{\link[=prob_statistic]{prob_statistic()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"} using \code{LF_1}.
#'
#' @param LF_2 Vector of statistics produced with
#'        \code{\link[=max_cong]{max_cong()}} or
#'        \code{\link[=max_incong]{max_incong()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"}.
#'
#' @param M02 Matrix produced with
#'        \code{\link[=prob_statistic]{prob_statistic()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"} using \code{LF_2}.
#'
#' @param LF_3 Vector of statistics produced with
#'        \code{\link[=max_cong]{max_cong()}} or
#'        \code{\link[=max_incong]{max_incong()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"}.
#'
#' @param M03 Matrix produced with
#'        \code{\link[=prob_statistic]{prob_statistic()}} for \code{"geoD"},
#'        \code{"paco"} or \code{"paraF"} using \code{LF_3}.
#'
#' @param ... Any optional argument admissible in
#'        \code{\link[graphics:boxplot]{boxplot()}}
#'
#' @return The Gini values obtained and their representation in a boxplot.
#'
#' @importFrom GiniWegNeg Gini_RSV
#' @importFrom graphics boxplot text
#' @export
#'
#' @examples
#' data(nuc_cp)
#' N = 10
#' n = 8
#' # Maximizing incongruence
#' NPi <- max_incong(np_matrix, NUCtr, CPtr, n, N, method = "paco",
#'                   symmetric = FALSE, ei.correct = "sqrt.D",
#'                   percentile = 0.99, diff.fq = TRUE,
#'                   strat = "parallel", cl = 8)
#' THSi <- trimHS_maxI(N, np_matrix, n)
#' PACOc <- prob_statistic(ths = THSi, np_matrix, NUC_500tr[1:5],
#'                         CP_500tr[1:5], freqfun = "paco", NPi,
#'                         symmetric = FALSE, ei.correct = "sqrt.D",
#'                         percentile = 0.99, diff.fq = TRUE, res.fq = FALSE,
#'                         below.p = TRUE, strat = "parallel", cl = 8)
gini_plot <- function (LF_1, M01, LF_2, M02, LF_3, M03,
                       ylab = "Normalized Gini coefficient", ...)
{
  if (missing(M01) & missing(LF_1) == TRUE) {
    LF01 <- NULL
    Gini01 <- NULL
    GiniM01 <- NULL
  }
  else {
    LF01 <- LF_1[, ncol(LF_1)]
    Gini01 <- unlist(Gini_RSV(LF_1))
    GiniM01 <- unlist(apply(M01, 1, Gini_RSV))
  }
  if (missing(M02) & missing(LF_2) == TRUE) {
    LF_2 <- NULL
    Gini02 <- NULL
    GiniM02 <- NULL
  }
  else {
    LF02 <- LF_2[, ncol(LF_2)]
    Gini02 <- unlist(Gini_RSV(LF_2))
    GiniM02 <- unlist(apply(M02, 1, Gini_RSV))
  }
  if (missing(M03) & missing(LF_3) == TRUE) {
    LF03 <- NULL
    Gini03 <- NULL
    GiniM03 <- NULL
  }
  else {
    LF03 <- LF_3[, ncol(LF_3)]
    Gini03 <- unlist(Gini_RSV(LF_3))
    GiniM03 <- unlist(apply(M03, 1, Gini_RSV))
  }
  ginis <- list(GiniM01, GiniM02, GiniM03)
  ginis <- ginis[!sapply(ginis, is.null)]
  rg <- list(Gini01 = Gini01, Gini02 = Gini02, Gini03 = Gini03)
  rg <- rg[!sapply(rg, is.null)]
  boxplot(ginis,
          show.names = TRUE, ...)
  for (i in 1:length(rg)) {
    text(i, rg[[i]], "*", cex = 2, col = "red")
  }
  return(unlist(rg))
}

