#' Selection frequency of each column across repeated selections
#'
#' @param sets A list of character vectors, one selected set per training
#'   sample (fold, resample or replicate).
#' @param columns All candidate column (or variable) names.
#' @return Named numeric vector: the share of `sets` that contain each column.
#' @export
selection_frequency <- function(sets, columns) {
  S <- vapply(sets, function(s) as.numeric(columns %in% s), numeric(length(columns)))
  stats::setNames(rowMeans(matrix(S, nrow = length(columns))), columns)
}

#' Stability index of Nogueira, Sechidis and Brown (2018)
#'
#' Summarizes selection frequencies over `M` selections among `p` columns.
#' It equals 1 when every selection gives the same set and has expected value
#' 0 when sets of the same average size are drawn at random.
#'
#' @param freq Selection frequency of each column, e.g. from
#'   [selection_frequency()]. Columns never selected may be left out if `p`
#'   is given.
#' @param M Number of selections.
#' @param p Total number of columns.
#' @return The stability index, or `NA` if it is undefined (fewer than two
#'   selections, or every selection empty or full).
#' @export
nogueira_index <- function(freq, M, p = length(freq)) {
  kbar <- sum(freq)
  if (M < 2 || kbar <= 0 || kbar >= p) return(NA_real_)
  s2 <- M / (M - 1) * freq * (1 - freq)
  1 - mean(c(s2, rep(0, p - length(freq)))) / ((kbar / p) * (1 - kbar / p))
}
