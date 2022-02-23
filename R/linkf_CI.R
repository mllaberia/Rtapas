#' Confidence intervals for the frequency of host-symbiont association
#'
#' From the matrix obtained in \code{\link[=prob_statistic]{prob_statistic()}},
#' compute the confidence intervals for the frequencies of the host-symbiont
#' associations using sets of posterior probability trees.
#'
#' @param freqfun Options are \code{"geo_D"}, \code{"paco_ss"} or
#'        \code{"paraF"}, depending on which confidence intervals you want
#'        to compute (apply to the result of \code{\link[=geo_D]{geo_D()}},
#'        \code{\link[=paco_ss]{paco_ss()}} or \code{\link[=paraF]{paraF()}}).
#'
#' @param x Matrix produced with \code{\link[=prob_statistic]{prob_statistic()}}
#'        for Geodesic distance, PACo or ParaFit.
#'
#' @param fx Vector of statistics produced with
#'        \code{\link[=link_freq]{link_freq()}} for Geodesic distance, PACo or
#'        ParaFit.
#'
#' @param c.level Confidence interval level. Default is \code{0.01} (99 per cent).
#'
#' @param barplot Default is \code{"TRUE"}, plots the distribution and confidence
#'        intervals of the frequencies.
#'
#' @param col.bar A vector of colors for the bars or bar components.
#'        By default, \code{"lightblue"} is used.
#'
#' @param col.ci A vector of colors for the confidence intervals arrows.
#'        By default, \code{"darkblue"} is used.
#'
#' @param ... Any graphical option admissible in
#'        \code{\link[=barplot]{barplot()}}
#'
#' @return A dataframe with the associations (columns 1 and 2), the observed
#'         value of the frequencies for these associations (column 3), the mean,
#'         the minimum and the maximum value of the frequencies (columns 4, 5
#'         and 6) obtained with the sets of posterior probability trees.
#'
#' @importFrom graphics arrows axis
#'
#' @export
#'
#' @examples
#' # ficus_aganoidae
#'
linkf_CI <- function (freqfun = "geo_D", x, fx, c.level = 95, barplot = TRUE,
                      col.bar = "lightblue", col.ci = "darkblue", ...) {
  freqfun.choice <- c("geoD", "paco", "paraF")
  if (freqfun %in% freqfun.choice == FALSE)
    stop(writeLines("Invalid freqfun parameter.\r Correct choices are 'geoD',\n
                    'paco' or 'paraF'"))
  if (freqfun == "geoD") {
    GD01 <- x
    LFGD01 <- fx
    a <- 1 - (c.level/100)
    GD.LO <- apply(GD01, 2, quantile, a/2)
    GD.HI <- apply(GD01, 2, quantile, 1 - (a/2))
    GD.AV <- apply(GD01, 2, mean)
    df <- data.frame(LFGD01[, 1], LFGD01[, 2], LFGD01[, ncol(LFGD01)],
                     GD.LO, GD.HI, GD.AV)
    colnames(df) <- c("Taxa1", "Taxa2", "GD.wFq",
                      "GD.LO", "GD.HI", "GD.AV")
    if (barplot == TRUE) {
      link.fq <- barplot(GD.AV, xaxt = "n", horiz = FALSE,
                         cex.names = 0.6, las = 2, cex.axis = 0.8,
                         ylab = "Observed - Expected frequency",
                         ylim = c(min(GD.LO), max(GD.HI)), col = col.bar,
                         ...)
      suppressWarnings(arrows(link.fq, GD.HI, link.fq,
                              GD.LO, length = 0, angle = 90, code = 3, col = col.ci))
      axis(side = 1, at = link.fq[1:length(GD.AV)], labels = LFGD01$HS,
           las = 2, tick = FALSE, line = 0.1, cex.axis = 0.5)
      return(df)
    }
    else {
      return(df)
    }
  }
  if (freqfun == "paco") {
    PACO01 <- x
    LFPACO01 <- fx
    a <- 1 - (c.level/100)
    PACO.LO <- apply(PACO01, 2, quantile, a/2)
    PACO.HI <- apply(PACO01, 2, quantile, 1 - (a/2))
    PACO.AV <- apply(PACO01, 2, mean)
    df <- data.frame(LFPACO01[, 1], LFPACO01[, 2], LFPACO01[, ncol(LFPACO01)], PACO.LO, PACO.HI, PACO.AV)
    colnames(df) <- c("Taxa1", "Taxa2", "PACO.wFq",
                      "PACO.LO", "PACO.HI", "PACO.AV")
    if (barplot == TRUE) {
      link.fq <- barplot(PACO.AV, xaxt = "n", horiz = FALSE,
                         cex.names = 0.6, las = 2, cex.axis = 0.8,
                         ylab = "Observed - Expected frequency",
                         ylim = c(min(PACO.LO), max(PACO.HI)), col = col.bar,
                         ...)
      suppressWarnings(arrows(link.fq, PACO.HI, link.fq,
                              PACO.LO, length = 0, angle = 90, code = 3, col = col.ci))
      axis(side = 1, at = link.fq[1:length(PACO.AV)], labels = LFPACO01$HS,
           las = 2, tick = FALSE, line = 0.1, cex.axis = 0.5)
      return(df)
    }
    else {
      return(df)
    }
  }
  if (freqfun == "paraF") {
    PF01 <- x
    LFPF01 <- fx
    a <- 1 - (c.level/100)
    PF.LO <- apply(PF01, 2, quantile, a/2)
    PF.HI <- apply(PF01, 2, quantile, 1 - (a/2))
    PF.AV <- apply(PF01, 2, mean)
    df <- data.frame(LFPF01[, 1], LFPF01[, 2], LFPF01[, ncol(LFPF01)],
                     PF.LO, PF.HI, PF.AV)
    colnames(df) <- c("Taxa1", "Taxa2", "PF.wFq",
                      "PF.LO", "PF.HI", "PF.AV")
    if (barplot == TRUE) {
      link.fq <- barplot(PF.AV, xaxt = "n", horiz = FALSE,
                         cex.names = 0.6, las = 2, cex.axis = 0.8,
                         ylab = "Observed - Expected frequency",
                         ylim = c(min(PF.LO), max(PF.HI)), col = col.bar,
                         ...)
      suppressWarnings(arrows(link.fq, PF.HI, link.fq,
                              PF.LO, length = 0, angle = 90, code = 3, col = col.ci))
      axis(side = 1, at = link.fq[1:length(PF.AV)], labels = LFPF01$HS,
           las = 2, tick = FALSE, line = 0.1, cex.axis = 0.5)
      return(df)
    }
    else {
      return(df)
    }
  }
}
