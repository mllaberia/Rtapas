#' Tanglegram of the host-symbiont frequencies
#'
#' Wrapper of \code{\link[phytools:plot.cophylo]{phytools::plot.cophylo()}}
#' is used for mapping as heatmap the host-symbiont frequencies estimated by
#' Random TaPas on a tanglegram. It also plots the average frequency
#' (or residual frequency) of occurrence of each terminal and optionally, the
#' fast maximum likelihood estimators of ancestral states of each node.
#'
#' @param treeH Host phylogeny. An object of class \code{"phylo"}.
#'
#' @param treeS Symbiont phylogeny. An object of class \code{"phylo"}.
#'
#' @param HS Host-symbiont association matrix.
#'
#' @param fqtab Dataframe produced with \code{\link[=link_freq]{link_freq()}}.
#'
#' @param colscale Choose between \code{"diverging"}, color reflects distance
#'        from 0 (centered at 0, recommended if \code{"res.fq = TRUE"})
#'        or \code{"sequential"}, color reflects distance from minimum value
#'        (spanning from the min to max frequencies observed).
#'
#' @param colgrad If \code{colscale =} vector defining the color gradient of the
#'        heatmap.
#'
#' @param nbreaks Number of discrete values along \code{"colgrad"}.
#'
#' @param node.tag Specifies whether maximum likelihood estimators of ancestral
#'        states are to be computed. Default is \code{TRUE}.
#'
#' @param cexpt Size of color points at terminals and nodes.
#'
#' @param link.lwd The line width for plotting, default to 1.
#'
#' @param link.lty The line type. Line types can either be specified as an
#'        integer (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted,
#'        4 = dotdash, 5 = longdash, 6 = twodash) or as one of the character
#'        strings \code{"blank"}, \code{"solid"}, \code{"dashed"},
#'        \code{"dotted"}, \code{"dotdash"}, \code{"longdash"}, or
#'        \code{"twodash"}. Note that \code{"blank"} uses invisible lines
#'        (i.e., does not draw them).
#'
#' @param fsize Relative font size for tip labels.
#'
#' @param pts Logical value indicating whether or not to plot filled circles at
#'        each vertex of the tree, as well as at transition points between
#'        mapped states. Default is \code{FALSE}.
#'
#' @param link.type If curved linking lines are desired, set to \code{"curved"}.
#'        Default is \code{"straight"}.
#'
#' @param ftype Font type. Options are \code{"reg"}, \code{"i"}
#'        (italics), \code{"b"} (bold) or \code{"bi"} (bold-italics).
#'
#' @param ... Any graphical option admissible in
#'        \code{\link[phytools:plot.cophylo]{phytools::plot.cophylo()}}
#'
#'
#' @return A tanglegram with quantitative information displayed as heatmap.
#'
#'
#' @examples
#' data(birds_mites)
#' N = 10
#' n = 50
#' TBM <- trimHS_maxC(N, bm_matrix, n, strat = "parallel", cl = 4)
#' PACO <- paco_ss(TBM, birds, mites, symmetric = TRUE,
#'                   proc.warns = FALSE, ei.correct = "sqrt.D",
#'                   strat = "parallel", cl = 8)
#' LFPACO <- link_freq(TBM, PACO, bm_matrix, percentile = 0.01, sep = "-",
#'                   below.p = TRUE, res.fq = TRUE)
#'
#' col.scale = c("darkred","gray90", "darkblue")
#' tangle_gram(birds, mites, bm_matrix, LFPACO, colscale= "diverging",
#'             colgrad = col.scale, nbreaks = 50, node.tag = FALSE,
#'             link.lty = 1, fsize = 0.25, pts = FALSE,
#'             link.type = "straight", cexpt = 1)
#'
#'
#' @import stats
#' @import phytools
#' @importFrom grDevices colorRampPalette
#'
#' @export
tangle_gram <- function(treeH, treeS, HS, fqtab, colscale = "diverging",
                        colgrad, nbreaks = 50, node.tag = TRUE, cexpt = 1,
                        link.lwd = 1, link.lty = 1, fsize = 0.5, pts = FALSE,
                        link.type = "straight", ftype = "i", ...) {

  colscale.choice <- c("diverging", "sequential")
  if (colscale %in% colscale.choice == FALSE)
    stop(writeLines("Invalid colscale parameter.\r
                     Correct choices are 'diverging' and 'sequential'"))

  colscale.range <- function(x) {
    rescale.range <- function(x) {
      xsq <- round(x)
      if(colscale == "sequential") {
        y <- range(xsq)
        col_lim <- (y[1]:y[2])-y[1]+1
        xsq <- xsq-y[1]+1
        new.range <- list(col_lim, xsq)
      } else {
        x1 <- x[which(x < 0)]
        if(length(x[which(x < 0)]) == 0) {
          warning("Not enough negative values for diverging scale.
                   The color scale is sequential")
          y <- range(xsq)
          col_lim <- (y[1]:y[2]) - y[1] + 1
          xsq <- xsq - y[1] + 1
          new.range <- list(col_lim, xsq)
        } else {
          x2 <- x[which(x >= 0)]
          x1 <- round(x1)
          x2 <- round(x2)
          y <- max(abs(x))
          col_lim <- (-y:y) + y + 1
          y1 <- range(x1)
          y2 <- range(x2)
          x1 <- x1-y1[1] + 1
          x2 <- x2-y2[1] + 1
          new.range <- list(col_lim, x1, x2)
        }
      }
      return(new.range)
    }
    if(colscale == "sequential" | length(x[which(x < 0)]) == 0) {
      NR <- rescale.range(x)
      rbPal <- colorRampPalette(colgrad)
      linkcolor <- rbPal(nbreaks)[as.numeric(cut(NR[[1]], breaks = nbreaks))]
      NR <- NR[[2]]
      linkcolor <- linkcolor[NR]
    } else {
      NR <- rescale.range(x)
      NR.neg <- NR[[1]] [which (NR[[1]] <= max(NR[[2]]))]
      NR.pos <- NR[[1]] [-NR.neg] - max(NR[[2]])
      m <- median(1:length(colgrad))
      colgrad.neg <- colgrad[which(1:length(colgrad) <= m)]
      colgrad.pos <- colgrad[which(1:length(colgrad) >= m)]
      rbPal <- colorRampPalette(colgrad.neg)
      linkcolor1 <- rbPal(nbreaks)[as.numeric(cut(NR.neg, breaks = nbreaks))]
      rbPal <- colorRampPalette(colgrad.pos)
      linkcolor2 <- rbPal(nbreaks)[as.numeric(cut(NR.pos, breaks = nbreaks))]
      linkcolor1 <- linkcolor1[NR[[2]]]
      linkcolor2 <- linkcolor2[NR[[3]]]
      linkcolor <- rep(NA, length(x))
      linkcolor[which(x < 0)] <- linkcolor1
      linkcolor[which(x >= 0)] <- linkcolor2
    }
    return(linkcolor)
  }

  if(colscale == "sequential"){
    links <- fqtab[, 4]
  } else {
    if(ncol(fqtab) < 5) {
      stop("'diverging' color scale requires residual frequencies
           (res.fq = TRUE) in link_f() function")
    } else {links <- fqtab[, 5]}
  }

  LKcolor <- colscale.range(links)
  HS.lut <- which(HS == 1, arr.ind = TRUE)
  linkhs <- cbind(rownames(HS)[HS.lut[, 1]], colnames(HS)[HS.lut[, 2]])
  obj <- phytools::cophylo(treeH, treeS, linkhs)
  phytools::plot.cophylo(obj, link.col = LKcolor, link.lwd = link.lwd,
                         link.lty = link.lty, fsize = fsize,
                         link.type = link.type, ftype = ftype, ...)

  Hfreq <- aggregate(links, by=list(freq = fqtab[, 1]), FUN=mean)
  Sfreq <- aggregate(links, by=list(freq = fqtab[, 2]), FUN=mean)

  Hfreq <- Hfreq[match(obj$trees[[1]]$tip.label, Hfreq$freq),]
  Sfreq <- Sfreq[match(obj$trees[[2]]$tip.label, Sfreq$freq),]

  if (node.tag == TRUE){
    fit.H <- fastAnc(obj$trees[[1]],Hfreq[, 2])
    fit.S <- fastAnc(obj$trees[[2]],Sfreq[, 2])
    NLH <- colscale.range (fit.H)
    NLS <- colscale.range (fit.S)
    nodelabels.cophylo(pch = 16, col = NLH, cex = cexpt)
    nodelabels.cophylo(pch = 16, col = NLS, cex = cexpt, which = "right")
  }
  TLH <- colscale.range (Hfreq[, 2])
  TLS <- colscale.range (Sfreq[, 2])
  phytools::tiplabels.cophylo(pch = 16, col = TLH, cex = cexpt)
  phytools::tiplabels.cophylo(pch = 16, col = TLS, cex = cexpt, which = "right")
}
