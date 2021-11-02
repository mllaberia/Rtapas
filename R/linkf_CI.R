#' Confidence intervals for the frequency of host-symbiont association
#'
#' @param ths List of trimmed matrices produced by \code{\link[=trimHS_maxC]{trimHS_maxC()}}.
#'
#' @param HS Host-Symbiont association matrix.
#'
#' @param mTreeH Number X of posterior-probabilistic trees of host.
#'
#' @param mTreeS Number X of posterior-probabilistic trees of symbiont.
#'
#' @param freqfun Options are \code{"geo_D"} or \code{"paco_ss"},
#'        depending on which confidence intervals you want to compute
#'        (apply to the result of \code{\link[=paco_ss]{paco_ss()}}
#'        or \code{\link[=geo_D]{geo_D()}}).
#'
#' @param geoD Vector of statistics produced with \code{\link[=geo_D]{geo_D()}}.
#'
#' @param paco Vector of statistics produced with \code{\link[=paco_ss]{paco_ss()}}.
#'
#' @param session Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"multisession"} and \code{"multicore"} (not
#'        supported on Windows) resolves \R expressions in parallel in separate
#'        \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with \code{\link[future:availableCores]{future::availableCores()}}.
#'        Note that \code{cl <= \link[=availableCores]{availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#' @param percentile Percentile to evaluate. Default is \code{0.01}. The percentile
#'        applied can be specified with \code{percentile.res.fq}, determines whether
#'        to apply a correction to the estimated frequencies by setting a null
#'        model in wich the occurrence of each host-symbiont association is
#'        evenly distributed along the whole frequency distribution.
#'
#' @param res.fq Determines whether a correction to avoid one-to-one associations
#'        being overrepresented in the percentile evaluated. If \code{TRUE} (default)
#'        a residual frequency value (observed - expected frequency) is computed
#'        for each host-symbiont association.
#'
#' @param below.p Determines whether frequencies are to be computed below or
#'        above the percentile set. Default is \code{TRUE}.
#'
#' @param symmetric Specifies the type of Procrustes superimposition. Default
#'        is \code{FALSE}, indicates that the superposition is applied
#'        asymmetrically (S depends on H). If \code{TRUE}, PACo is applied
#'        symmetrically (dependency between S and H is reciprocal).
#'
#' @param ei.correct Specifies how to correct potential negative eigenvalues
#'        from the conversion of phylogenetic distances into Principal
#'        Coordinates: \code{"none"} (the default) indicates that no correction
#'        is required, particularly if H and S are ultrametric; \code{"sqrt.D"}
#'        takes the element-wise square-root of the phylogenetic distances;
#'        \code{"lingoes"} and \code{"cailliez"} apply the classical Lingoes and
#'        Cailliez corrections, respectively.
#'
#' @param barplot Default is \code{"TRUE"}, plots the distribution and confidence
#'        intervals of the frequencies.
#'
#' @param ... Any graphical option admissible in \code{\link[=barplot]{barplot()}}

#'
#' @return A dataframe with
#'
#' @export
#'
#' @examples
linkf_CI <- function (ths, HS, mTreeH, mTreeS, freqfun = "geo_D", geoD, paco,
                      percentile = 0.01, res.fq = TRUE, below.p = TRUE, symmetric=FALSE, ei.correct="none",
                      session = "sequential", cl = 1, barplot = TRUE, ...) {

  freqfun.choice <- c("geo_D", "paco_ss")
  if(freqfun %in% freqfun.choice == FALSE)
  stop(writeLines("Invalid freqfun parameter.\r Correct choices are 'geo_D' or 'paco_ss'"))

  if(freqfun == "geo_D") {
    LFGD01 <- link_freq(ths, geoD, HS, percentile = percentile,
                        res.fq = res.fq, below.p = below.p)
    GD01 <- matrix(NA, length(mTreeH), nrow(LFGD01))
    if(session == "sequential"){
      for(i in 1:length(mTreeH)) {
        GD.CI <- geo_D(ths, treeH=mTreeH[[i]], treeS=mTreeS[[i]], session = session, cl = cl)
        LFGD01.CI <- link_freq(ths, GD.CI, HS, percentile = percentile,
                               res.fq = res.fq, below.p = below.p)
        GD01[i,] <- LFGD01.CI[,5]
        }
      } else {
        for(i in 1:length(mTreeH)) {
          GD.CI <- geo_D(ths, treeH=mTreeH[[i]], treeS=mTreeS[[i]], session = session, cl = cl)
          LFGD01.CI <- link_freq(ths, GD.CI, HS, percentile = percentile,
                                 res.fq = res.fq, below.p = below.p)
          GD01[i,] <- LFGD01.CI[,5]
          }
      }
    colnames(GD01) <- LFGD01[,3]
    GD.LO <- apply(GD01, 2, quantile, 0.025)
    GD.HI <- apply(GD01, 2, quantile, 0.975)
    GD.AV <- apply(GD01, 2, mean)

    df <- data.frame(GD01[i,], GD.LO, GD.HI, GD.AV)
    colnames(df) <- c("GDwFq", "GD.LO", "GD.HI", "GD.AV")

    if(barplot == TRUE) {
      link.fq <-barplot(GD.AV, xaxt='n',
                        horiz=FALSE, cex.names = 0.6, las=2, cex.axis=0.8,
                        ylab="Observed - Expected frequency",
                        ylim=c(min(GD.LO), max(GD.HI)), col="lightblue")
      suppressWarnings(arrows(link.fq, GD.HI, link.fq, GD.LO, length= 0,
                              angle=90, code=3, col="darkblue"))
      axis(side=1, at=link.fq[1:length(GD.AV)], labels=LFGD01$HS, las=2,
           tick = FALSE, line= 0.1, cex.axis=0.5)

      return(df)
    } else {return(df)}
  }

  if(freqfun == "paco_ss") {
    LFPACO01 <- link_freq (ths, paco, HS, percentile = percentile,
                           res.fq = res.fq, below.p = below.p)
    PACO01 <- matrix(NA, length(mTreeH), nrow(LFPACO01))
    if(session == "sequential"){
      for(i in 1:length(mTreeH)) {
        PA.CI <- paco_ss(ths, treeH=mTreeH[[i]], treeS=mTreeS[[i]],
                         symmetric = symmetric, ei.correct = ei.correct,
                         session = session, cl = cl)
        LFPA01.CI <- link_freq(ths, PA.CI, HS, percentile = percentile,
                               res.fq = res.fq, below.p = below.p)
        PACO01[i,] <- LFPA01.CI[,5]
        }
      } else {
        for(i in 1:length(mTreeH)) {
          PA.CI <- paco_ss(ths, treeH=mTreeH[[i]], treeS=mTreeS[[i]],
                           symmetric = symmetric, ei.correct = ei.correct,
                           session = session, cl = cl)
          LFPA01.CI <- link_freq(ths, PA.CI, HS, percentile = percentile,
                                 res.fq = res.fq, below.p = below.p)
          PACO01[i,] <- LFPA01.CI[,5]
          }
      }
    colnames(PACO01) <- LFPACO01[,3]
    PACO.LO <- apply(PACO01, 2, quantile, 0.025)
    PACO.HI <- apply(PACO01, 2, quantile, 0.975)
    PACO.AV <- apply(PACO01, 2, mean)

    if(barplot == TRUE){
      link.fq <-barplot(PACO.AV, xaxt='n',
                        horiz=FALSE, cex.names = 0.6, las=2, cex.axis=0.8,
                        ylab="Observed - Expected frequency",
                        ylim=c(min(PACO.LO), max(PACO.HI)),col="lightblue")
      suppressWarnings(arrows(link.fq, PACO.HI, link.fq, PACO.LO, length= 0,
                              angle=90, code=3, col="darkblue"))
      axis(side=1, at=link.fq[1:length(PACO.AV)], labels=LFPACO01$HS, las=2,
           tick = FALSE, line= 5, cex.axis=0.6)

      return(PACO.AV)
    } else {return(PACO.AV)}
  }
}
