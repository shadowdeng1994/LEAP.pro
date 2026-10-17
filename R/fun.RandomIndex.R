#' Generate a short random index for output file names (internal).
#'
#' Drawn from the global RNG stream, so `set.seed()` makes it reproducible.
#'
#' @param n Number of characters. (default: 4)
#'
#' @return A length-1 character string of `n` characters from `[0-9A-Za-z]`.
#' @noRd
fun.RandomIndex <- function(n=4){
  tmp.pool <- c(0:9,letters,LETTERS)
  paste0(sample(tmp.pool,n,replace=TRUE),collapse="")
}
