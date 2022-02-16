#' Plot the confidence intervals of Gini coefficient
#'
#' Computes and plot the normalized Gini coefficient (G*) (Raffinetti et al.
#' 2015) and its confidence intervals of the residual frequency distributions
#' of the  desired statistics: \code{"geo_D"}, \code{"paco_ss"} or
#' \code{"paraF"} (either all or only one of them).
#'
#' @param geoD Matrix produced with \code{\link[=prob_statistic]{prob_statistic()}}
#'        for Geodesic distance.
#' @param lkgd Vector of statistics produced with
#'        \code{\link[=link_freq]{link_freq()}} for Geodesic distance.
#' @param paco Matrix produced with \code{\link[=prob_statistic]{prob_statistic()}}
#'        for PACo.
#' @param lkpaco Vector of statistics produced with
#'        \code{\link[=link_freq]{link_freq()}} for PACo.
#' @param parafit Matrix produced with \code{\link[=prob_statistic]{prob_statistic()}}
#'        for ParaFit.
#' @param lkpf Vector of statistics produced with
#'        \code{\link[=link_freq]{link_freq()}} for ParaFit.
#' @param ...
#'
#' @return
#'
#' @importFrom GiniWegNeg Gini_RSV
#' @export
#'
#' @examples
gini_plot <- function(geoD, lkgd, paco, lkpaco, parafit, lkpf, ...) {

  if(missing(geoD) & missing(lkgd) == TRUE) {
    lkgd <- NULL
    GiniGD <- NULL
    GiniMGD <- NULL
  } else {
    LFGD <- lkgd[, ncol(lkgd)]
    GiniGD <- unlist(Gini_RSV(LFGD))
    GiniMGD <- unlist(apply(geoD, 1, Gini_RSV))
  }
  if(missing(paco) & missing(lkpaco) == TRUE) {
    lkpaco <- NULL
    GiniPA <- NULL
    GiniMPA <- NULL
  } else {
    LFPACO <- lkpaco[, ncol(lkpaco)]
    GiniPA <- unlist(Gini_RSV(LFPACO))
    GiniMPA <- unlist(apply(paco, 1, Gini_RSV))
  }
  if(missing(parafit) & missing(lkpf) == TRUE) {
    lkpf <- NULL
    GiniPF <- NULL
    GiniMPF <- NULL
  } else {
    LFPF <- lkpf[, ncol(lkpf)]
    GiniPF <- unlist(Gini_RSV(LFPF))
    GiniMPF <- unlist(apply(parafit, 1, Gini_RSV))
  }

  ginis <- list(GiniMGD, GiniMPA, GiniMPF)
  ginis <- ginis[!sapply(ginis, is.null)]

  rg <- list("G_GD" = GiniGD, "G_PA" = GiniPA, "G_PF" = GiniPF)
  #rg <- list(GiniGD, GiniPA, GiniPF)
  rg <- rg[!sapply(rg, is.null)]

  boxplot(ginis, ylab="Normalized Gini coefficient",
          las=3, show.names = TRUE, ...)
      for (i in 1:length(rg)) {
        text(i, rg[[i]],"*", cex = 2, col = "red")
      }
  return(unlist(rg))
}

