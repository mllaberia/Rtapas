#' Plot the confidence intervals of Gini coefficient
#'
#' Computes and displays in a boxplot the Gini coefficient and their
#' confidence intervals of the frequency (or residual/corrected frequencies)
#' distributions of up to three estimated (in)congruence metrics of the
#' individual host-symbiont associations.
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
#' @param ylab Title of the y label.
#'
#' @return The Gini values obtained and their representation in a boxplot, with
#'         their confidence intervals.
#'
#' @section NOTE:
#'          It produces a conventional Gini coefficient (G)
#'          (Ultsch and Lötsch 2017) if all output values are positive, or
#'          a normalized Gini coefficient (G*) (Raffinetti et al. 2015) if
#'          negative values are produced due to corrected frequencies
#'          (if \code{res.fq = TRUE} or
#'          \code{diff.fq = TRUE}). For more details see
#'          \code{\link[GiniWegNeg:Gini_RSV]{Gini_RSV()}}.
#'
#' @references
#' Ultsch A., Lötsch J. (2017). A data science based standardized Gini index
#' as a Lorenz dominance preserving measure of the inequality of distributions.
#' PLOS ONE. 12:e0181572. \doi{10.1371/journal.pone.0181572}
#'
#' Raffinetti E., Siletti E., Vernizzi A. (2015). On the Gini coefficient
#' normalization when attributes with negative values are considered.
#' Stat Methods Appl. 24:507–521. \doi{10.1007/s10260-014-0293-4}
#'
#'
#'
#'
#' @importFrom GiniWegNeg Gini_RSV
#' @importFrom graphics boxplot text
#' @export
#'
#' @examples
#' data(nuc_cp)
#' N = 10000
#' n = 15
#' # Maximizing congruence
#' NPc_PACo <- max_cong(np_matrix, NUCtr, CPtr, n, N, method = "paco",
#'                 symmetric = FALSE, ei.correct = "sqrt.D",
#'                 percentile = 0.01, res.fq = FALSE,
#'                 strat = "parallel", cl = 8)
#'
#' # Loaded directly from dataset
#' # THSC <- trimHS_maxC(N, np_matrix, n)
#' # pp_treesPACo_cong <- prob_statistic(ths = THSc, np_matrix, NUC_500tr[1:10],
#' #                         CP_500tr[1:10], freqfun = "paco", NPc,
#' #                         symmetric = FALSE, ei.correct = "sqrt.D",
#' #                         percentile = 0.01, res.fq = FALSE,  below.p = TRUE,
#' #                         strat = "parallel", cl = 8)
#'
#' gini_plot(LF_1 = NPc_PACo, M01 = pp_treesPACo_cong,
#'           ylab = "Gini Coefficient (G)",
#'           names = c("PACo_MaxCongruence"))
#'
gini_plot <- function (LF_1, M01, LF_2, M02, LF_3, M03,
                       ylab = "Gini coefficient", ...)
{
  if (missing(M01) & missing(LF_1) == TRUE) {
    Gini01 <- NULL
    GiniM01 <- NULL
  }
  else {
    Gini01 <- unlist(Gini_RSV( LF_1[, ncol(LF_1)]))
    GiniM01 <- unlist(apply(M01, 1, Gini_RSV))
  }
  if (missing(M02) & missing(LF_2) == TRUE) {
    Gini02 <- NULL
    GiniM02 <- NULL
  }
  else {
    Gini02 <- unlist(Gini_RSV(LF_2[, ncol(LF_2)]))
    GiniM02 <- unlist(apply(M02, 1, Gini_RSV))
  }
  if (missing(M03) & missing(LF_3) == TRUE) {
    Gini03 <- NULL
    GiniM03 <- NULL
  }
  else {
    Gini03 <- unlist(Gini_RSV(LF_3[, ncol(LF_3)]))
    GiniM03 <- unlist(apply(M03, 1, Gini_RSV))
  }
  ginis <- list(GiniM01, GiniM02, GiniM03)
  ginis <- ginis[!sapply(ginis, is.null)]
  rg <- list(Gini01 = Gini01, Gini02 = Gini02, Gini03 = Gini03)
  rg <- rg[!sapply(rg, is.null)]
  boxplot(ginis,
          show.names = TRUE, ylab = ylab, ...)
  for (i in 1:length(rg)) {
    text(i, rg[[i]], "*", cex = 2, col = "red")
  }
  return(unlist(rg))
}

