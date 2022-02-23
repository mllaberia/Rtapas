#' Algortihm for maximizing congruence between two phylogenies
#'
#' @param HS Host-Symbiont association matrix.
#'
#' @param treeH Host phyolgeny. An object of class "phylo".
#'
#' @param treeS Symbiont phylogeny. An object of class "phylo".
#'
#' @param n Number of unique associations.
#'
#' @param N Number of runs.
#'
#' @param method cc
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
#'        \code{"lingoes"} and \code{"cailliez"} apply the classical Lingoes
#'        and Cailliez corrections, respectively.
#'
#' @param percentile Percentile to evaluate. Default is \code{0.01}.
#'
#' @param res.fq Determines whether a correction to avoid one-to-one
#'        associations being overrepresented in the percentile evaluated.
#'        If \code{TRUE} (default) a residual frequency value (observed -
#'        expected frequency) is computed for each host-symbiont association.
#'
#' @param diff.fq FALSE
#'
#' @param strat Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"parallel"} resolves \R expressions in parallel in
#'        separate \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with
#'        \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#'
#' @return A dataframe with host-symbiont associations in rows. The first and
#'         second columns display the names of the host and symbiont terminals,
#'         respectively. The third column designates the host-symbiont
#'         association by pasting the names of the terminals, and the fourth
#'         column displays the frequency of occurrence of each host-symbiont
#'         association. If \code{res.fq = TRUE}, column 5 displays the
#'         corrected frequencies as a residual
#'
#' @export
#'
#'
#' @examples
#' diff.fq
#'
max_cong <- function (HS, treeH, treeS, n, N, method = "paco", symmetric = FALSE,
                        ei.correct = "none", percentile = 0.01, res.fq = TRUE,
                        diff.fq = FALSE, strat = "sequential", cl = 1){

  THS <- trimHS_maxC(N = N, n = n, HS = HS, check.unique = TRUE)
  method.choice <- c("geoD", "paco", "paraF")
  if (method %in% method.choice == FALSE)
    stop(writeLines("Invalid global-fit method. Correct choices are 'paco' or 'parafit'"))
  if (method == "geoD") {
    GD <- geo_D(THS, treeH, treeS, strat = strat, cl = cl)

    LF <- link_freq(THS, GD, HS, percentile = percentile, below.p = TRUE,
                    res.fq = res.fq)

    if (diff.fq == TRUE) {
      LFi_b <- link_freq(THS, GD, HS, percentile = 0.01,
                         below.p = TRUE, res.fq = res.fq)
      LFi_u <- link_freq(THS, GD, HS, percentile = 0.99,
                         below.p = FALSE, res.fq = res.fq)
      LFr <- LFi_b[, ncol(LFi_b)] - LFi_u[, ncol(LFi_u)]
      LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
      return(LFi_n)
    } else {return(LF)}
  }
  if (method == "paco") {
    PACO <- paco_ss(THS, treeH, treeS, symmetric = symmetric,
                    proc.warns = FALSE, ei.correct = ei.correct,
                    strat = strat, cl = cl)
    LF <- link_freq(THS, PACO, HS, percentile = percentile, below.p = TRUE,
                    res.fq = res.fq)

    if (diff.fq == TRUE) {
      LFi_b <- link_freq(THS, PACO, HS, percentile = 0.01,
                         below.p = TRUE, res.fq = res.fq)
      LFi_u <- link_freq(THS, PACO, HS, percentile = 0.99,
                         below.p = FALSE, res.fq = res.fq)
      LFr <- LFi_b[, ncol(LFi_b)] - LFi_u[, ncol(LFi_u)]
      LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
      return(LFi_n)
    } else {return(LF)}
  }

  if (method == "paraF") {
    PF <- paraF(THS, treeH, treeS, ei.correct = ei.correct,
                strat = strat, cl = cl)
    LF <- link_freq(THS, PF, HS, percentile = percentile, below.p = TRUE,
                    res.fq = res.fq)

    if (diff.fq == TRUE) {
      LFi_b <- link_freq(THS, PF, HS, percentile = 0.01,
                         below.p = TRUE, res.fq = res.fq)
      LFi_u <- link_freq(THS, PF, HS, percentile = 0.99,
                         below.p = FALSE, res.fq = res.fq)
      LFr <- LFi_b[, ncol(LFi_b)] - LFi_u[, ncol(LFi_u)]
      LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
      return(LFi_n)
    } else {return(LF)}
  }
}
