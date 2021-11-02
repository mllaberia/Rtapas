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
#' @param session Strategy you want to work with. Default is \code{"sequential"},
#'        resolves \R expressions sequentially in the current \R
#'        process. If \code{"multisession"} and \code{"multicore"} (not
#'        supported on Windows) resolves \R expressions in parallel in separate
#'        \R sessions running in the background.
#'
#' @param cl Number of cluster the user wants to use. Check how many CPUs/cores
#'        your computer has with \code{\link[future:availableCores]{future::availableCores()}}.
#'        Note that \code{cl<=\link[=availableCores]{availableCores()}}.
#'        Default is \code{cl = 1} for \code{"sequential"} strategy.
#'
#'
#' @return A list of the N trimmed matrices.
#'
#' @examples
#' #trimHS_maxC(N, HS, n, check.unique = TRUE)
#'
#' @importFrom future plan availableCores
#' @import future.apply
#'
#' @export
#'
trimHS_maxC <- function (N, HS, n, check.unique = TRUE,
                           session = "sequential", cl = 1) {

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

  session.choice <- c("sequential", "multisession", "multicore")
  if (session %in% session.choice == FALSE)
    stop(writeLines("Invalid session parameter"))

  if (session == "sequential") {
    plan(sequential)
    trim.HS <- future_lapply(1:N, trim.int, HS = HS, n = n, future.seed = NULL)
    if (check.unique == TRUE) trim.HS <- unique(trim.HS)
    if (length(trim.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trim.HS[sapply(trim.HS, is.null)] <- NULL
    return(trim.HS)
  } else {
    plan(session, workers = cl)
    trim.HS <- future_lapply(1:N, trim.int, HS = HS, n = n, future.seed = NULL)
    plan(sequential)
    if (check.unique == TRUE) trim.HS <- unique(trim.HS)
    if (length(trim.HS) < N)
      warning("No. of trimmed H-S assoc. matrices < No. of runs")
    trim.HS[sapply(trim.HS, is.null)] <- NULL
    return(trim.HS)
  }
}


