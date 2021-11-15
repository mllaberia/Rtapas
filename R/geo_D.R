#' Geodesic distance between trees
#'
#' For any trimmed matrix produced with \code{\link[=trimHS_maxC]{trimHS_maxC()}},
#' it prunes the host-symbiont phylogenies to conform with the trimmed matrix
#' and computes the geodesic distance between the pruned trees.
#' \code{NOTE}: This function can only be used with strictly bifurcating trees.
#'
#' @param ths Trimmed matrix.
#'
#' @param treeH Host phyolgeny. An object of class \code{"phylo"}.
#'
#' @param treeS Symbiont phylogeny. An object of class \code{"phylo"}.
#'
#' @param strat Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"parallel"} resolves \R expressions in parallel in
#'        separate \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#' @return Geodesic distance
#'
#' @examples
#' # birds_mites dataset
#' data(birds_mites)
#' N = 1e+2
#' n = 50
#' TBM <- trimHS_maxC(N, bm_matrix, n, strat = "parallel", cl = 4)
#' gd_bm <- geo_D(TBM, treeH = birds, treeS = mites, strat = "parallel", cl = 8)
#'
#' # plant_fungi dataset
#' data(plant_fungi)
#' N = 1e+2
#' n = 15
#' TPF <- trimHS_maxC(N, pf_matrix, n)
#' gd_pf <- geo_D(TPF, treeH = plant, treeS = fungi, strat = "parallel", cl = 8)
#'
#'
#' @import ape
#' @import distory
#'
#' @export
geo_D <- function(ths, treeH, treeS,
                  strat = "sequential", cl = 1) {

  geoD <- function (ths, treeH, treeS) {
    treeh <- drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
    trees <- drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))

    ths <- ths[treeh$tip.label, trees$tip.label]

    ths.lut <- which(ths[treeh$tip.label, trees$tip.label]==1, arr.ind = TRUE)
    dummy.labels <- rownames(ths.lut)
    trees$tip.label <- dummy.labels
    combo.tree <- list(treeh, trees)
    gd <- dist.multiPhylo(combo.tree)
    return(gd)
  }

  strat.choice <- c("sequential", "parallel")
  if (strat %in% strat.choice == FALSE)
    stop(writeLines("Invalid strategy parameter"))

  if (strat == "sequential") {
    gd_f <- sapply(ths, geoD, treeH = treeH, treeS = treeS)
    return(gd_f)
  } else {
    cores <- makeClusterPSOCK(workers = cl)
    gd_f <- parSapply(cores, ths, geoD, treeH = treeH, treeS = treeS)
    stopCluster(cores)
    return(gd_f)
  }
}
