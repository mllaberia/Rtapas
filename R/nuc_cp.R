#' @docType data
#'
#' @title Nuclear and chloroplast dataset of orchids
#' @name nuc_cp
#' @usage data(nuc_cp)
#'
#' @description Data set of nuclear and chloroplast loci of 52 orchid taxa
#'              from Kew DNA and Tissue Collection,
#'              https://dnabank.science.kew.org/homepage.html
#'              (Perez-Escobar et al. 2021).
#'
#' @format This data set compresses five objects:
#'
#' \describe{
#'    \item{\code{np_matrix}}{Associations one-to-one between the 52 orchid
#'    taxa.
#'    A binary matrix with 52 rows (nuclear) and 52 columns (chloroplast).}
#'
#'    \item{\code{NUCtr}}{Phylogeny constructed by sequence data of nuclear
#'    loci of orchids (Perez-Escobar et al. 2021).
#'    An object of class \code{"phylo"} containing a list with the details of
#'    the phylogenetic tree (i.e. edge, edge length, nodes and tips names).}
#'
#'    \item{\code{CPtr}}{Phylogeny constructed by sequence data of chloroplast
#'    loci of orchids (Perez-Escobar et al. 2021).
#'    An object of class \code{"phylo"} containing a list with the details of
#'    the phylogenetic tree (i.e. edge, edge length, nodes and tips names).}
#'
#'    \item{\code{NUC_500tr}}{500 bootstrap replicates ultrametric trees
#'    from Perez-Escobar et al. (2021).
#'    List of class \code{"multiphylo"} containing a 500 phylogenetic trees
#'    with their respective details (i.e. edges, edges length, nodes, and
#'    tips names).}
#'
#'    \item{\code{CP_500tr}}{500 bootstrap replicates ultrametric trees
#'    from Perez-Escobar et al. (2021).
#'    List of class \code{"multiphylo"} containing a 500 phylogenetic trees
#'    with their respective details (i.e. edges, edges length, nodes, and
#'    tips names).}
#'  }
#'
#' @references
#' Perez-Escobar O.A., Dodsworth S., Bogarin D., Bellot S., Balbuena J.A.,
#' Schley R., Kikuchi I., Morris S.K., Epitawalage N., Cowan R., Maurin O.,
#' Zuntini A., Arias T., Serna A., Gravendeel B., Torres M.F., Nargar K.,
#' Chomicki G., Chase M.W., Leitch I.J., Forest F., Baker W.J. (2021).
#' Hundreds of nuclear and plastid loci yield novel insights into orchid
#' relationships. American Journal of Botany, 108(7), 1166-1180.
#' \doi{10.1002/ajb2.1702}
#'
#' @source
#' Perez-Escobar O.A., Dodsworth S., Bogarin D., Bellot S., Balbuena J.A.,
#' Schley R., Kikuchi I., Morris S.K., Epitawalage N., Cowan R., Maurin O.,
#' Zuntini A., Arias T., Serna A., Gravendeel B., Torres M.F., Nargar K.,
#' Chomicki G., Chase M.W., Leitch I.J., Forest F., Baker W.J. (2021).
#' Hundreds of nuclear and plastid loci yield novel insights into orchid
#' relationships. American Journal of Botany, 108(7), 1166-1180.
#' \doi{10.1002/ajb2.1702}
#'
#' @keywords datasets
#'
#'
NULL
