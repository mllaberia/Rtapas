#' Create an host-symbiont association matrix
#'
#' @param hs
#'
#' @return
#' @export
#'
#' @examples
ass_m <- function(hs) {
  host <- hs[,1]
  symb <- hs[,2]
  mhs <- (matrix(NA, nrow = length(symb), ncol = length(host)))
  diag(mhs) <- 1
  mhs[is.na(mhs)] <- 0
  colnames(mhs) <- symb
  rownames(mhs) <- host

  m <- mhs %>% by(rownames(mhs), colSums)
  m1 <- as.matrix(do.call(rbind,m))

  m2 <- t(m1)

  m2 <- m2 %>% by(rownames(m2), colSums)
  m2 <- as.matrix(do.call(cbind,m2))
  return(m2)
}
