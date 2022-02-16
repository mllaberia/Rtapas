#' @docType data
#'
#' @title birds_mites dataset
#' @name birds_mites
#' @usage data(birds_mites)
#'
#' @description Data set of 93 passerine birds species and and their feather
#' mites.
#'
#' @format This data set compresses four objects:
#'
#' \describe{
#'    \item{\code{bm_matrix}}{Associations between 93 passerine birds
#'    species and their feather mite, collected in eight regions
#'    (Costa Rica, Kazakhstan, Mexico, Panama, Peru, Russia, Tanzania,
#'    USA) by Klimov et al. (2017).
#'    A binary matrix with 93 rows (birds) and 98 variables (feather mites).}
#'
#'    \item{\code{birds}}{Phylogeny of passerine brids extracted from
#'    one of 200 random stationary Bayesian time-calibrated trees
#'    downloaded form "A global phylogeny of birds" http://birdtree.org.
#'    An object of class \code{"phylo"} containing a list with the details of
#'    the phylogenetic tree (i.e. edges, edges length, nodes, root edge and
#'    tips names).}
#'
#'    \item{\code{mites}}{Consensus tree of feather mites constructed by
#'    Klimov et al. (2017).
#'    An object of class \code{"phylo"} containing a list with the details of
#'    the phylogenetic tree (i.e. edges, edges length, nodes and tips names).}
#'
#'    \item{\code{birds_200r_trees}}{200 random stationary Bayesian
#'    time-calibrated trees downloaded form "A global phylogeny of birds"
#'    http://birdtree.org.
#'    List of class \code{"multiphylo"} containing a 200 phylogenetic trees
#'    with their respective details (i.e. edges, edges length, nodes, root edge
#'    and tips names).}
#'  }
#'
#' @references
#' Klimov, P.B., Mironov, S.V., OConnor, B.M. (2017). Detecting ancient
#' codispersals and host shifts by double dating of host and parasite
#' phylogenies: Application in proctophyllodid feather mites associated with
#' passerine birds. Evolution, 71(10), 2381-2397. \doi{10.1111/evo.13309}
#'
#' @source
#' Balbuena, J. A., Pérez-Escobar, A.A., Llopis-Belenguer, C.,
#' Blasco-Costa, I. (2020). Random tanglegram partitions (Random TaPas):
#' An alexandrian approach to the cophylogenetic gordian knot. Systematic
#' biology, 69(6), 1212-1230. \doi{10.1093/sysbio/syaa033}.
#'
#' @keywords datasets
#'
#'
NULL

