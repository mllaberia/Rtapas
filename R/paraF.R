#' Test of host-symbiont coevolution
#'
#' @param ths Trimmed matrix.
#'
#' @param treeH Host phylogeny. An object of class \code{"phylo"}.
#'
#' @param treeS Symbiont phylogeny. An object of class \code{"phylo"}.
#'
#' @param ei.correct Specifies how to correct potential negative eigenvalues
#'        from the conversion of phylogenetic distances into Principal
#'        Coordinates: \code{"none"} (the default) indicates that no correction
#'        is required, particularly if H and S are ultrametric; \code{"sqrt.D"}
#'        takes the element-wise square-root of the phylogenetic distances;
#'        \code{"lingoes"} and \code{"cailliez"} apply the classical Lingoes and
#'        Cailliez corrections, respectively.
#'
#' @param strat Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"parallel"} resolves \R expressions in parallel in
#'        separate \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy
#'
#' @return A number object with the statistic of the global host-symbiont test
#'         for any trimmed matrix.
#'
#' @import ape
#' @export
#'
#' @examples
#' #paraF()
paraF <- function (ths, treeH, treeS, ei.correct = "none",
                     strat = "sequential", cl = 1) {
  paraf <- function (ths, treeH, treeS, ...) {
    eigen.choice <- c("none", "lingoes", "cailliez", "sqrt.D")
    if (ei.correct %in% eigen.choice == FALSE)
      stop(writeLines("Invalid eigenvalue correction parameter.\r
               Correct choices are 'none', 'lingoes', 'cailliez' or 'sqrt.D'"))
    treeh <- ape::drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
    trees <- ape::drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))
    if (is.null(treeh)==TRUE | is.null(trees)==TRUE) PF <- NA else {
    # Reorder ths as per tree labels:
    ths <- ths[treeh$tip.label, trees$tip.label]
    DH <- ape::cophenetic.phylo(treeh)
    DP <- ape::cophenetic.phylo(trees)
    if (ei.correct == "sqrt.D"){DH <- sqrt(DH); DP <- sqrt(DP); ei.correct ="none"}
    PF <- ape::parafit(DH, DP, ths, nperm=1, silent=TRUE)
    PF <- PF$ParaFitGlobal
    }
    return(PF)
  }

  strat.choice <- c("sequential", "parallel")
  if (strat %in% strat.choice == FALSE)
    stop(writeLines("Invalid strategy parameter"))

  if (strat == "sequential") {
    pfit <- sapply(ths, paraf, treeH, treeS, ei.correct = ei.correct)
    return(pfit)
  } else {
    cores <- makeClusterPSOCK(workers = cl)
    pfit <- parSapply(cores, ths, paraf, treeH, treeS, ei.correct = ei.correct)
    stopCluster(cores)
    return(pfit)
  }
}
