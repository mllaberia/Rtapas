#' Test of host-symbiont coevolution
#'
#' For any trimmed matrix produced with
#' \code{\link[=trimHS_maxC]{trimHS_maxC()}} or
#' \code{\link[=trimHS_maxI]{trimHS_maxI()}}, it prunes the host (H) & symbiont
#' (S) phylogenies to conform with the trimmed matrix and runs
#' \code{\link[ape:parafit]{ape::parafit()}} (Legendre et al. 2002) to
#' calculate the ParaFitGlobal Statistic.
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
#'        your computer has with
#'        \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy
#'
#' @return A number object with the ParaFitGlobal Statistic of host-symbiont
#'         test for any trimmed matrix.
#'
#' @import ape
#' @import parallel
#' @importFrom parallelly makeClusterPSOCK
#'
#'
#' @export
#'
#' @examples
#' # birds_mites dataset
#' data(birds_mites)
#' N = 1e+2
#' n = 50
#' TBM <- trimHS_maxC(N, bm_matrix, n, strat = "parallel", cl = 4)
#' PARAF <- pacoF(TBM, birds, mites, ei.correct = "sqrt.D",
#'                   strat = "parallel", cl = 8)
#'
paraF <- function (ths, treeH, treeS, ei.correct = "none",
                     strat = "sequential", cl = 1) {
  paraf <- function (ths, treeH, treeS, ...) {
    eigen.choice <- c("none", "lingoes", "cailliez", "sqrt.D")
    if (ei.correct %in% eigen.choice == FALSE)
      stop(writeLines("Invalid eigenvalue correction parameter.\r
               Correct choices are 'none', 'lingoes', 'cailliez' or 'sqrt.D'"))
    treeh <- drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
    trees <- drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))
    if (is.null(treeh)==TRUE | is.null(trees)==TRUE) PF <- NA else {
    # Reorder ths as per tree labels:
    ths <- ths[treeh$tip.label, trees$tip.label]
    DH <- cophenetic.phylo(treeh)
    DP <- cophenetic.phylo(trees)
    if (ei.correct == "sqrt.D"){DH <- sqrt(DH); DP <- sqrt(DP); ei.correct ="none"}
    PF <- parafit(DH, DP, ths, nperm = 0, silent = TRUE, correction = ei.correct)
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
    a <- length(ths)/cl
    b <- ceiling(a)
    ths1 <- split(ths, rep(1:cl, each = b))
    pf <- c()
    cores <- makeClusterPSOCK(workers = cl)
    for(i in 1:length(ths1)) {
      ths2 <- ths1[[i]]
      pfit <- parSapply(cores, ths2, paraf, treeH, treeS, ei.correct = "none")
      pf <- c(pf, pfit)
    }
    stopCluster(cores)
    return(pf)
  }
}
