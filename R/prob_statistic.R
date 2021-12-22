#' Frequencies of the associations for the posterior probabilistic trees
#'
#' Computes the frequencies of the associations for each of the probabilistic
#' trees from the statistics of \code{\link[=geo_D]{geo_D()}},
#' \code{\link[=paco_ss]{paco_ss()}} or \code{\link[=paraF]{paraF()}}).
#'
#' @param ths List of trimmed matrices produced by
#'        \code{\link[=trimHS_maxC]{trimHS_maxC()}}.
#'
#' @param HS Host-Symbiont association matrix.
#'
#' @param mTreeH Number X of posterior-probabilistic trees of host.
#'
#' @param mTreeS Number X of posterior-probabilistic trees of symbiont.
#'
#' @param freqfun Options are \code{"geo_D"}, \code{"paco_ss"} or
#'        \code{"paraF"}, depending on which confidence intervals you want to
#'        compute (apply to the result of \code{\link[=geo_D]{geo_D()}},
#'        \code{\link[=paco_ss]{paco_ss()}} or \code{\link[=paraF]{paraF()}})
#'
#' @param fx Vector of statistics produced with
#'         \code{\link[=link_freq]{link_freq()}} for Geodesic distance, PACo or
#'         ParaFit.
#'
#' @param percentile Percentile to evaluate. Default is \code{0.01}.
#'
#' @param below.p Determines whether frequencies are to be computed below or
#'        above the percentile set. Default is \code{TRUE}.
#'
#' @param res.fq Determines whether a correction to avoid one-to-one
#'        associations being overrepresented in the percentile evaluated.
#'        If \code{TRUE} (default) a residual frequency value (observed -
#'        expected frequency) is computed for each host-symbiont association.

#' @param symmetric Specifies the type of Procrustes superimposition. Default
#'        is \code{FALSE}, indicates that the superposition is applied
#'        asymmetrically (S depends on H). If \code{TRUE}, PACo is applied
#'        symmetrically (dependency between S and H is reciprocal).
#'
#' @param proc.warns Switches on/off trivial warnings returned when treeH and
#'        treeS differ in size (number of tips). Default is \code{FALSE}.
#'
#' @param ei.correct Specifies how to correct potential negative eigenvalues
#'        from the conversion of phylogenetic distances into Principal
#'        Coordinates: \code{"none"} (the default) indicates that no correction
#'        is required, particularly if H and S are ultrametric; \code{"sqrt.D"}
#'        takes the element-wise square-root of the phylogenetic distances;
#'        \code{"lingoes"} and \code{"cailliez"} apply the classical Lingoes
#'        and Cailliez corrections, respectively.
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
#' @return A matrix with the value of the statistics for each of the
#'         probabilistic trees.
#' @export
#'
#' @examples
prob_statistic <- function (ths, HS, mTreeH, mTreeS, freqfun = "geo_D", fx,
                       percentile = 0.01, res.fq = TRUE, below.p = TRUE,
                       symmetric=FALSE, ei.correct="none", proc.warns = FALSE,
                       strat = "sequential", cl = 1) {

  strat.choice <- c("sequential", "parallel")
  if (strat %in% strat.choice == FALSE)
    stop(writeLines("Invalid strategy parameter"))

  freqfun.choice <- c("geo_D", "paco_ss", "paraF")
  if(freqfun %in% freqfun.choice == FALSE)
    stop(writeLines("Invalid freqfun parameter.\r Correct choices are 'geo_D',
                    'paco_ss' or 'paraF'"))

  if(freqfun == "geo_D") {
    LFGD01 <- link_freq(ths, fx, HS, percentile = percentile,
                        res.fq = res.fq, below.p = below.p)

    gd_apply <- function(ths, mTreeH, mTreeS, HS) {
      GD.CI <- geo_D(ths, treeH=mTreeH, treeS=mTreeS)
      LFGD01.CI <- link_freq(ths, GD.CI, HS, percentile = percentile,
                             res.fq = res.fq, below.p = below.p)
      return(LFGD01.CI[, ncol(LFGD01.CI)])
    }
    if(strat == "sequential") {
    GD01 <- t(sapply(1:length(mTreeH), function(x) gd_apply(ths, mTreeH[[x]], mTreeS[[x]], HS)))
    colnames(GD01) <- LFGD01[,3]
    return(GD01)
    } else {
      geoD <- function (ths, treeH, treeS) {
        treeh <- ape::drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
        trees <- ape::drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))

        ths <- ths[treeh$tip.label, trees$tip.label]

        ths.lut <- which(ths[treeh$tip.label, trees$tip.label]==1, arr.ind = TRUE)
        dummy.labels <- rownames(ths.lut)
        trees$tip.label <- dummy.labels
        combo.tree <- list(treeh, trees)
        gd <- distory::dist.multiPhylo(combo.tree)
        return(gd)
      }
      cores <- makeClusterPSOCK(workers = cl)
      GD01 <- matrix(NA, length(mTreeH), nrow(LFGD01))
      for(i in 1:length(mTreeH)) {
        GD.CI <- parallel::parSapply(cores, ths, geoD, treeH=mTreeH[[i]],
                                     treeS= mTreeS[[i]])
        LFGD01.CI <- RandomTaPas::link_freq(ths, GD.CI, HS, percentile = percentile,
                                            res.fq = res.fq, below.p = below.p)
        GD01[i,] <- LFGD01.CI[,5]
      }
      stopCluster(cores)
      colnames(GD01) <- LFGD01[,3]
    }
    return(GD01)
    } else {
      if(freqfun == "paco_ss") {
        LFPACO01 <- link_freq (ths, fx, HS, percentile = percentile,
                             res.fq = res.fq, below.p = below.p)

        paco_apply <- function(ths, mTreeH, mTreeS, HS) {
          PA.CI <- paco_ss(ths, treeH=mTreeH, treeS=mTreeS,
                         symmetric = symmetric, ei.correct = ei.correct,
                         strat = strat , cl = cl)
          LFPA01.CI <- link_freq(ths, PA.CI, HS, percentile = percentile,
                               res.fq = res.fq, below.p = below.p)
          return(LFPA01.CI[, ncol(LFPA01.CI)])
        }
        if(strat == "sequential") {
          PACO01 <- t(sapply(1:length(mTreeH), function(x) paco_apply(ths, mTreeH[[x]], mTreeS[[x]], HS)))
          colnames(PACO01) <- LFPACO01[,3]
          return(PACO01)
          } else {
            pacoss <- function (ths, treeH, treeS, ...) {
          eigen.choice <- c("none", "lingoes", "cailliez", "sqrt.D")
          if (ei.correct %in% eigen.choice == FALSE)
            stop(writeLines("Invalid eigenvalue correction parameter.\r
               Correct choices are 'none', 'lingoes', 'cailliez' or 'sqrt.D'"))
          treeh <- ape::drop.tip(treeH, setdiff(treeH$tip.label, rownames(ths)))
          trees <- ape::drop.tip(treeS, setdiff(treeS$tip.label, colnames(ths)))
          # Reorder ths as per tree labels:
          ths <- ths[treeh$tip.label, trees$tip.label]
          DH <- ape::cophenetic.phylo(treeh)
          DP <- ape::cophenetic.phylo(trees)
          if(ei.correct == "sqrt.D"){DH <- sqrt(DH); DP <- sqrt(DP); ei.correct ="none"}
          D <- paco::prepare_paco_data(DH, DP, ths)
          D <- paco::add_pcoord(D, correction = ei.correct)
          if (proc.warns == FALSE) D <- vegan::procrustes(D$H_PCo, D$P_PCo,
                                                        symmetric = symmetric)
          else
            D <- suppressWarnings(vegan::procrustes(D$H_PCo, D$P_PCo,
                                                  symmetric = symmetric))
          return(D$ss)
          }
            cores <- makeClusterPSOCK(workers = cl)
            PACO01 <- matrix(NA, length(mTreeH), nrow(LFPACO01))
            for(i in 1:length(mTreeH)) {
              PA.CI<-parallel::parSapply(cores, ths, pacoss, treeH=mTreeH[[i]],
                                   treeS= mTreeS[[i]], symmetric = symmetric,
                                   proc.warns = proc.warns, ei.correct = ei.correct)
              LFPA01.CI <- link_freq(ths, PA.CI, HS,  percentile = percentile,
                               res.fq = res.fq, below.p = below.p)
              PACO01[i,] <- LFPA01.CI[,5]
            }
            stopCluster(cores)
            colnames(PACO01) <- LFPACO01[,3]
            }
        return(PACO01)

      } else {
        LFPF01 <- link_freq (ths, fx, HS, percentile = percentile,
                               res.fq = res.fq, below.p = below.p)
        pf_apply <- function(ths, mTreeH, mTreeS, HS) {
          PF.CI <- paraF(ths, treeH=mTreeH, treeS=mTreeS,
                           ei.correct = ei.correct,
                           strat = strat , cl = cl)
          LFPF01.CI <- link_freq(ths, PF.CI, HS, percentile = percentile,
                                 res.fq = res.fq, below.p = below.p)
          return(LFPF01.CI[, ncol(LFPF01.CI)])
        }
        if(strat == "sequential") {
          PF01 <- t(sapply(1:length(mTreeH), function(x) pf_apply(ths, mTreeH[[x]], mTreeS[[x]], HS)))
          colnames(PF01) <- LFPF01[,3]
          return(PF01)
        } else {
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
              PF <- ape::parafit(DH, DP, ths, nperm=1, silent=TRUE, correction = ei.correct)
              PF <- PF$ParaFitGlobal
            }
            return(PF)
          }
          cores <- makeClusterPSOCK(workers = cl)
          PF01 <- matrix(NA, length(mTreeH), nrow(LFPF01))
          for(i in 1:length(mTreeH)) {
            PF.CI<-parallel::parSapply(cores, ths, paraf, treeH=mTreeH[[i]],
                                       treeS= mTreeS[[i]],
                                       proc.warns = proc.warns, ei.correct = ei.correct)
            LFPF01.CI <- link_freq(ths, PF.CI, HS,  percentile = percentile,
                                   res.fq = res.fq, below.p = below.p)
            PF01[i,] <- LFPF01.CI[,5]
          }
          stopCluster(cores)
          colnames(PF01) <- LFPF01[,3]
        }
        return(PF01)


      }
    }
}

