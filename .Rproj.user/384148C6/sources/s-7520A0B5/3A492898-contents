#' Frequency of host-symbiont association
#'
#' After applying \code{\link[=geo_D]{geo_D()}} or \code{\link[=paco_ss]{paco_ss()}}
#' to each matrix produced by \code{\link[=trimHS_maxC]{trimHS_maxC()}}, this
#' function determines the frequency of each host-symbiont association occurring
#' in a given percentile of cases that maximize phylogenetic congruence (either
#' the fraction of lowest geodesic distances or lowest sum of squared residuals).
#'
#' @param x List of trimmed matrices produced by \code{\link[=trimHS_maxC]{trimHS_maxC()}}.
#'
#' @param fx Vector of statistics produced with either
#'         \code{\link[=geo_D]{geo_D()}} or \code{\link[=paco_ss]{paco_ss()}}.
#'
#' @param HS Host-symbiont association matrix.
#'
#' @param percentile Percentile to evaluate. Default is \code{0.01}. The percentile
#'        applied can be specified with \code{percentile.res.fq}, determines whether
#'        to apply a correction to the estimated frequencies by setting a null
#'        model in wich the occurrence of each host-symbiont association is
#'        evenly distributed along the whole frequency distribution.
#'
#' @param sep Character that separates host and symbiont labels.
#'
#' @param below.p Determines whether frequencies are to be computed below or
#'        above the percentile set. Default is \code{TRUE}.
#'
#' @param res.fq Determines whether a correction to avoid one-to-one associations
#'        being overrepresented in the percentile evaluated. If \code{TRUE} (default)
#'        a residual frequency value (observed - expected frequency) is computed
#'        for each host-symbiont association.
#'
#' @section NOTE:
#'        The \code{res.fq = TRUE} correction is recommended in tanglegrams with
#'        large portion of multiple (as opposed to one-to-one) host-symbiont
#'        associations. For future usage, frequencies of host-symbiont
#'        associations above a given percentile values can also be computed
#'        setting \code{below.p = FALSE}.
#'
#'
#' @return A dataframe with host-symbiont associations in rows. The first and
#'         second columns display the names of the host and symbiont terminals,
#'         respectively. The third column designates the host-symbiont
#'         association by pasting the names of the terminals, and the fourth
#'         column displays the frequency of occurrence of each host-symbiont
#'         association. If \code{res.fq = TRUE}, column 5 displays the corrected
#'         frequencies as a residual.
#'
#' @examples
#' #link_freq
#'
#' @export
link_freq <- function (x, fx, HS, percentile = 0.01,
                       sep= "-", below.p = TRUE, res.fq = TRUE) {
  if (below.p == TRUE)
    percent <- which(fx <= quantile(fx, percentile, na.rm = TRUE)) else
      percent <- which(fx >= quantile(fx, percentile, na.rm = TRUE))
    trim.HS <- x[percent]
    paste.link.names <- function(X, sep) {
      X.bin <- which(X>0, arr.ind = TRUE)
      Y <- diag(nrow(X.bin))
      Y <- diag(nrow(X.bin))
      rownames(Y) <- rownames(X)[X.bin[,1]]
      colnames(Y) <- colnames(X)[X.bin[,2]]
      pln <- paste(rownames(Y), colnames(Y), sep = sep)
      return(pln)
    }
    link.names <- t(sapply(trim.HS, paste.link.names, sep = sep))
    lf <- as.data.frame(table(link.names))
    HS.LUT <- which(HS == 1, arr.ind = TRUE)
    linkf <- as.data.frame(cbind(rownames(HS)[HS.LUT[,1]],
                                 colnames(HS)[HS.LUT[,2]]))
    colnames(linkf) <- c('H', 'S')
    linkf$HS <- paste(linkf[,1], linkf[,2], sep = sep)
    linkf$Freq <- rep(0, nrow(linkf))
    linkf[match(lf[,1],linkf[,3]), 4] <- lf[,2]
    linkf2 <- linkf
    #
    if (res.fq == TRUE) {
      link.names.all <- t(sapply(x, paste.link.names, sep = sep))
      lf.all <- as.data.frame(table(link.names.all))
      linkf.all <- as.data.frame(cbind(rownames(HS)[HS.LUT[,1]],
                                       colnames(HS)[HS.LUT[,2]]))
      colnames(linkf.all) <- c('H', 'S')
      linkf.all$HS <- paste(linkf.all[,1], linkf.all[,2], sep = sep)
      linkf.all$Freq <- rep(0, nrow(linkf.all))
      linkf.all[match(lf.all[,1], linkf.all[,3]), 4] <- lf.all[,2]
      w <- linkf.all[,4]
      w <- as.matrix(w*percentile)
      wFq <- linkf$Freq - w
      linkf$wFq <- wFq
    } else linkf <- linkf2
    return(linkf)
}
