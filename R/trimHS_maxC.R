#' Trims the H-S association matrix maximizing the Congruence
#'
#' For N runs, it randomly chooses \code{n} unique one-to-one associations and
#' trims the H-S association matrix to include only the n associations.
#'
#' @param N Number of runs.
#'
#' @param HS Host-Symbiont association matrix.
#'
#' @param n Number of unique associations.
#'
#' @param check.unique if \code{TRUE} discards duplicated trimmed matrices. This
#'        alternative is recommended if \code{n} is small, because the
#'        probability of obtaining the same trimmed matrix in different runs
#'        increases as \code{n} decreases.
#'
#' @param strat Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"parallel"} resolves \R expressions in parallel in
#'        separate \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with \code{\link[parallelly:availableCores]{parallelly::availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#'
#' @return A list of the N trimmed matrices.
#'
#' @examples
#' # With birds_mites dataset
#' data(birds_mites)
#' N = 1e+2
#' n = 50
#' TBM <- trimHS_maxC(N, bm_matrix, n, check.unique = TRUE, strat = "parallel", cl = 4)
#'
#' # With plant_fungi dataset
#' data(plant_fungi)
#' N = 1e+2
#' n = 15
#' TPF <- trimHS_maxC(N, pf_matrix, n)
#'
#'
#' @import parallelly
#' @import parallel
#'
#' @export
#'
trimHS_maxC <- function (N, HS, n, check.unique = TRUE,
                           strat = "sequential", cl = 1) {

  trim.int <- function (x, HS, n) {
    HS.LUT <- which(HS == 1, arr.ind = TRUE)
    HS.LUT <- cbind(HS.LUT, 1:nrow(HS.LUT))
    df <- as.data.frame(HS.LUT)
    hs.lut <- subset(df[sample(nrow(df)), ],
                     !duplicated(row) & !duplicated(col))
    if (nrow(hs.lut) < n) hs <- NULL else {
      hs.lut <- hs.lut[sample(nrow(hs.lut), n), ]
      hs <- diag(nrow(hs.lut))
      rownames(hs) <- rownames(HS[hs.lut[ ,1], ])
      colnames(hs) <- colnames(HS[ ,hs.lut[ ,2]])
      return(hs)
    }
  }

  strat.choice <- c("sequential", "parallel")
  if (strat %in% strat.choice == FALSE)
    stop(writeLines("Invalid strategy parameter"))

  if (strat == "sequential") {
    trim.HS <- lapply(1:N, trim.int, HS = HS, n = n)
    if (check.unique == TRUE) trim.HS <- unique(trim.HS)
    if (length(trim.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trim.HS[sapply(trim.HS, is.null)] <- NULL
    return(trim.HS)
  } else {
    cores <- makeClusterPSOCK(workers = cl, autoStop = TRUE)
    trim.HS <- parLapply(cores, 1:N, trim.int, HS = HS, n = n)
    if (check.unique == TRUE) trim.HS <- unique(trim.HS)
    if (length(trim.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trim.HS[sapply(trim.HS, is.null)] <- NULL
    return(trim.HS)
  }
}


