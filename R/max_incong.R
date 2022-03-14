#' Algortihm for maximizing incongruence between two phylogenies
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
#' @param method Specifies the desired global-fit method (GD, PACo or ParaFit).
#'        The default is \code{PACo}.
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
#' @param percentile Percentile to evaluate (\emph{p}). Default is
#'        \code{0.99} (99%).
#'
#' @param diff.fq Determines whether a correction to detect those associations
#'        that present higher incongruence in the system.
#'        If \code{TRUE} a corrected frequency value (observed in p -
#'        observed in (p-1)) is computed for each host-symbiont association.
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
#' @return A dataframe with host-symbiont associations in rows. The first and
#'         second columns display the names of the host and symbiont terminals,
#'         respectively. The third column designates the host-symbiont
#'         association by pasting the names of the terminals, and the fourth
#'         column displays the frequency of occurrence of each host-symbiont
#'         association in \emph{p}. If \code{diff.fq = TRUE}, column 5 displays
#'         the corrected frequencies.
#'
#' @section NOTE:
#'       If the \code{node.label} object in both trees contains NAs or neither
#'       value, an error will appear. You should make sure that all nodes
#'       have a value, or else remove it within the \code{"phylo"} class tree
#'       with \code{tree$node.label <- NULL}. For more details, see
#'       \code{\link[distory:dist.multiPhylo]{distory::dist.multiPhylo()}}
#'
#'       \code{GD} method can not be used with the trimmed matrices produced
#'       with \code{\link[=trimHS_maxI]{trimHS_maxI()}} or with the algorithm
#'       \code{\link[=max_incong]{max_incong()}} for those datasets with
#'       multiple associations.
#'
#' @export
#'
#' @examples
#' data(nuc_pc)
#' N = 1e+2
#' n = 15
#' NPi <- max_incong(np_matrix, NUCtr, CPtr, n, N, method = "paco",
#'                   symmetric = FALSE, ei.correct = "sqrt.D",
#'                   percentile = 0.99, diff.fq = TRUE,
#'                   strat = "parallel", cl = 10)
#'
#'
max_incong <- function (HS, treeH, treeS, n, N, method = "paco",
                        symmetric = FALSE, ei.correct = "none",
                        percentile = 0.99, diff.fq = FALSE,
                        strat = "sequential", cl = 1) {

  THSi <- trimHS_maxI(N = N, n = n, HS = HS, check.unique = TRUE)
  method.choice <- c("geoD", "paco", "paraF")
  if (method %in% method.choice == FALSE)
    stop(writeLines("Invalid global-fit method.
                    Correct choices are 'geoD', 'paco' or 'parafit'"))

  if (method == "geoD") {
    if(sum(HS) == ncol(HS)){
      GD <- geo_D(THSi, treeH, treeS, strat = strat, cl = cl)

      LF <- link_freq(THSi, GD, HS, percentile = percentile, below.p = FALSE,
                      res.fq = res.fq)

      if (diff.fq == TRUE) {
        LFi_b <- link_freq(THSi, GD, HS, percentile = 0.01,
                           below.p = TRUE, res.fq = FALSE)
        LFi_u <- link_freq(THSi, GD, HS, percentile = 0.99,
                           below.p = FALSE, res.fq = FALSE)
        LFr <- LFi_b[, ncol(LFi_b)] - LFi_u[, ncol(LFi_u)]
        LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
        return(LFi_n)
      } else {
        return(LF)
      }
      } else {stop("All associations must be one-to-one to apply the incongruent
               algorithm with this method due to the dissimilarity between the
               trimmed matrices")}
  }

  if (method == "paco") {
    PACO <- paco_ss(THSi, treeH, treeS, symmetric = symmetric,
                    proc.warns = FALSE, ei.correct = ei.correct, strat = strat,
                    cl = cl)
    LF <- link_freq(THSi, PACO, HS, percentile = percentile,
                    below.p = FALSE, res.fq = FALSE)

    if (diff.fq == TRUE) {
      LFi_b <- link_freq(THSi, PACO, HS, percentile = 0.01,
                         below.p = TRUE, res.fq = FALSE)
      LFi_u <- link_freq(THSi, PACO, HS, percentile = 0.99,
                         below.p = FALSE, res.fq = FALSE)
      LFr <- LFi_b[, ncol(LFi_b)] - LFi_u[, ncol(LFi_u)]
      LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
      return(LFi_n)
    } else {return(LF)}
  }

  if (method == "paraF") {
    PF <- paraF(THSi, treeH, treeS, ei.correct = ei.correct,
                strat = strat, cl = cl)
    LF <- link_freq(THSi, PF, HS, percentile = percentile,
                    below.p = FALSE, res.fq = FALSE)

    if (diff.fq == TRUE) {
      LFi_b <- link_freq(THSi, PF, HS, percentile = 0.01,
                         below.p = TRUE, res.fq = FALSE)
      LFi_u <- link_freq(THSi, PF, HS, percentile = 0.99,
                         below.p = FALSE, res.fq = FALSE)
      LFr <- LFi_u[, ncol(LFi_u)] - LFi_b[, ncol(LFi_b)]
      LFi_n <- cbind(LFi_b[, 1:ncol(LFi_b)], LFr)
      return(LFi_n)
    } else {return(LF)}
  }
}

