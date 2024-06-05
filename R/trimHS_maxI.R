#' Trims the H-S association matrix maximizing the incongruence
#'
#' For N runs, it randomly chooses \code{n} associations and trims the H-S
#' association matrix to include them, allowing both single and multiple
#' associations.
#'
#' @param N Number of runs.
#'
#' @param HS Host-Symbiont association matrix.
#'
#' @param n Number of associations.
#'
#' @param check.unique if \code{TRUE} discards duplicated trimmed matrices. This
#'        alternative is recommended if \code{n} is small, because the
#'        probability of obtaining the same trimmed matrix in different runs
#'        increases as \code{n} decreases.
#'
#' @param strat Flag indicating whether execution is to be  \code{"sequential"}
#'        or \code{"parallel"}. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"parallel"} resolves \R expressions in parallel in
#'        separate \R sessions running in the background.
#'
#' @param cl Number of cluster to be used for parallel computing.
#'        \code{\link[parallelly:availableCores]{parallelly::availableCores()}}
#'        returns the number of clusters available.
#'        Default is \code{cl = 1} resulting in \code{"sequential"} execution.
#'
#' @return A list of the N trimmed matrices.
#'
#' @export
#'
#' @examples
#' # data(nuc_cp)
#' # N = 10  #for the example, we recommend 1e+4 value
#' # n = 15
#' # TNC <- trimHS_maxI(N, np_matrix, n, check.unique = TRUE)
#'
trimHS_maxI <- function (N, HS, n, check.unique = TRUE,
                         strat = "sequential", cl = 1) {

  trim.intI <- function (x, HS, n) {
    HS.LUT <- which(HS == 1, arr.ind = TRUE)
    HS.LUT <- cbind(HS.LUT, 1:nrow(HS.LUT))
    LH <- LS <- 2
    while(LH <= 2 | LS <= 2) {  # skip configurations with less than
      HS.trim <- HS.LUT[sample(nrow(HS.LUT), n), ] # three H or three S
      LH <- length(unique(HS.trim[, 1]))
      LS <- length(unique(HS.trim[, 2]))
    }
    hs <- HS
    hs[ , ] <- 0
    for (j in 1:nrow(HS.trim)) {
      hs[HS.trim[j,1], HS.trim[j,2]] <- 1
    }
    hs <- hs[-which(rowSums(hs) == 0), -which(colSums(hs) == 0)]
    return(hs)
  }

  strat.choice <- c("sequential", "parallel")
  if (strat %in% strat.choice == FALSE)
    stop(writeLines("Invalid strategy parameter"))

  if (strat == "sequential") {
    trimI.HS <- lapply(1:N, trim.intI, HS = HS, n = n)
    if (check.unique == TRUE) trimI.HS <- unique(trimI.HS)
    if (length(trimI.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trimI.HS[sapply(trimI.HS, is.null)] <- NULL

  } else {
    cores <- parallelly::makeClusterPSOCK(workers = cl)
    trimI.HS <- parallel::parLapply(cores, 1:N, trim.intI, HS = HS, n = n)
    parallel::stopCluster(cores)
    if (check.unique == TRUE) trimI.HS <- unique(trimI.HS)
    if (length(trimI.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trimI.HS[sapply(trimI.HS, is.null)] <- NULL
  }
  dims <- function(x){all(dim(x)) == 0}
  which_null <- sapply(trimI.HS, dims)
  if (length(which(which_null == TRUE)) > 1) {
    trimI.HS <- trimI.HS[-which(which_null == TRUE)]
    }
  return(trimI.HS)
}

