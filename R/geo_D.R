#' Geodesic distance between trees
#'
#' For any trimmed matrix produced with
#' \code{\link[=trimHS_maxC]{trimHS_maxC()}} it prunes the host-symbiont
#' phylogenies to conform with the trimmed matrix and computes the geodesic
#' distance between the pruned trees.
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
#'        your computer has with
#'        \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#' @section NOTE:
#'       If the \code{node.label} object in both trees contains NAs or neither
#'       value, an error will appear. You should make sure that all nodes
#'       have a value, or else remove it within the \code{"phylo"} class tree
#'       with \code{tree$node.label <- NULL}. For more details, see
#'       \code{\link[distory:dist.multiPhylo]{distory::dist.multiPhylo()}}
#'
#'       This function can not be used with the trimmed matrices produced
#'       with \code{\link[=trimHS_maxI]{trimHS_maxI()}}.
#'
#' @return Geodesic distance
#'
#' @examples
#' # birds_mites dataset
#'
#' @import ape
#' @import distory
#' @import parallel
#' @importFrom parallelly makeClusterPSOCK
#'
#' @export
geo_D <- function(ths, treeH, treeS,
                  strat = "sequential", cl = 1) {

  if(is.binary(treeH) == FALSE | is.binary(treeS) == FALSE)
    stop("All trees must be strictly binary")

  geoD <- function (ths, treeH, treeS) {
    treeh <- drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
    trees <- drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))

    ths <- ths[treeh$tip.label, trees$tip.label]

    ths.lut <- which(ths[treeh$tip.label, trees$tip.label] == 1, arr.ind = TRUE)
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
