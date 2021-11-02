#' Maximum n for pick unique one-to-one association over a number of runs
#'
#' For a matrix of host-symbiont associations, it finds the maximum \code{n} for
#' which one-to-one unique associations can be picked in
#' \code{\link[=trimHS_maxC]{trimHS_maxC()}} over a number of runs.
#'
#' @param HS Host-symbiont association matrix.
#'
#' @param reps Number of runs to evaluate.
#'
#' @return The maximum number of unique associations (\code{n}).
#'
#' @section NOTE:
#'          It can be used to decide the best \code{n} prior to application of
#'          \code{\link[=trimHS_maxC]{trimHS_maxC()}}.
#'
#' @import stringr
#'
#' @examples
#' #one2one_f
#'
#' @export
one2one_f <- function(HS, reps = 1e+4, session = "sequential", cl = 1,
                      plot = TRUE, na.percent = 0.2){

  one2one <- function (HS, ...) {
    HS.LUT <- which(HS == 1, arr.ind = TRUE)
    HS.LUT <- cbind(HS.LUT,1:nrow(HS.LUT))
    df <- as.data.frame(HS.LUT)
    V <- rep(NA,reps)
    for(i in 1:reps){
      hs.lut <- subset(df[sample(nrow(df)),],
                       !duplicated(row) & !duplicated(col))
      n <- sum(HS)
      while (n >0) {
        n <- n-1;
        if (nrow(hs.lut) == n) break
      }
      V[i] <- n
    }
    V <- min(V)
    return(V)
  }
  np <- function(c, d){
    c <- na.percent
    d <- round(sum(HS) * c)
    return(d)
  }
  one  <- one2one(HS, reps)
  r.one <- round(one * 0.1)
  a <- (one - r.one):(one + r.one)
  if (plot == TRUE) {
    b <- rep(NA, length(a))
    for (i in 1:length(a)) {
      THS <- trimHS_maxC(reps, HS, n=a[i], check.unique = TRUE, session, cl)
      b[i] <- length(THS)
    }
    plot(a, b, type = "b", xlim = c(min(a), max(a)), xaxt = "n",
         xlab = "Number of unique H-S associations",
         ylab = "Number of runs accomplished") +
      axis(1, seq(round(min(a)), round(max(a)), by = 1),
           labels = (min(a):max(a)))

    x <- c("n.max", "n.rank", "posible.n")
    y <- c(one, stringr::str_c("[", 1, ";", one , "]"), np(c,d))
    r <- data.frame(x, y)
  } else {
    x <- c("n.max", "n.rank", "posible.n")
    y <- c(one, stringr::str_c("[", 1, ";", one , "]"), np(c,d))
    r <- data.frame(x, y)
  }
  return(r)
}

