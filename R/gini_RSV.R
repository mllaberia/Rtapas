#' The Gini coefficient adjusted for negative attributes (Raffinetti, Siletti,
#' & Vernizzi, 2015)
#'
#' Computes the Gini coefficient adjusted for negative (even weighted) data.
#'
#' @param y a vector of attributes containing even negative elements
#'
#' @section NOTE:
#'          It produces a conventional Gini coefficient (G)
#'          (Ultsch and Lötsch 2017) if all output values are positive, or
#'          a normalized Gini coefficient (G*) (Raffinetti et al. 2015) if
#'          negative values are produced due to corrected frequencies
#'          (if \code{res.fq = TRUE} or
#'          \code{diff.fq = TRUE}). For more details see
#'          Raffinetti et al. (2015).
#'
#' @references
#' Ultsch A., Lötsch J. (2017). A data science based standardized Gini index
#' as a Lorenz dominance preserving measure of the inequality of distributions.
#' PLOS ONE. 12:e0181572. \doi{10.1371/journal.pone.0181572}
#'
#' Raffinetti E., Siletti E., Vernizzi A. (2015). On the Gini coefficient
#' normalization when attributes with negative values are considered.
#' Stat Methods Appl. 24:507–521. \doi{10.1007/s10260-014-0293-4}
#'
#' @export
#'
#' @examples
#' data(nuc_cp)
#' N = 10 #for the example, we recommend 1e+4 value
#' n = 15
#' # Maximizing congruence
#' NPc_PACo <- max_cong(np_matrix, NUCtr, CPtr, n, N, method = "paco",
#'                 symmetric = FALSE, ei.correct = "sqrt.D",
#'                 percentile = 0.01, res.fq = FALSE,
#'                 strat = "parallel", cl = 8)
#'
#' gini_RSV(y = NPc_PACo)
#'
gini_RSV <- function(y){
  w <- rep(1, length(y))
  if (is.numeric(y) == TRUE){y = y} else {y <- y[, ncol(y)]}
  dataset<-cbind(y,w)
  ord_y<-order(y)
  dataset_ord<-dataset[ord_y,]
  y<-dataset_ord[,1]
  w<-dataset_ord[,2]
  N<-sum(w)
  yw<-y*w
  C_i<-cumsum(w)
  num_1<-sum(yw*C_i)
  num_2<-sum(yw)
  num_3<-sum(yw*w)
  G_num<-(2/N^2)*num_1-(1/N)*num_2-(1/N^2)*num_3
  t_neg<-subset(yw,yw<=0)
  T_neg<-sum(t_neg)
  T_pos<-sum(yw)+abs(T_neg)
  n_RSV<-(2*(T_pos+(abs(T_neg)))/N)
  mean_RSV<-(n_RSV/2)
  G_RSV<-(1/mean_RSV)*G_num
  GINI_RSV=G_RSV
  return(GINI_RSV)
}

