#' Creates a host-symbiont association matrix
#'
#' @param hs A two-columns matrix representing associations between host
#' (column 1) and symbiont (column 2) species.
#'
#' @return An association matrix, with hosts in rows and symbionts
#' in columns, sorted alphabetically.
#' @export
#'
#' @examples
ass_m <- function(hs) {
  host <- hs[,1]
  symb <- hs[,2]
  mhs <- (matrix(NA, nrow = length(host), ncol = length(symb)))
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
