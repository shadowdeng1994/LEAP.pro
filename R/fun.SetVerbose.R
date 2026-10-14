#' @title Set the LEAP log verbosity level.
#' @description
#' Controls the amount of progress output produced by the LEAP functions (e.g.
#' \code{GetInternalMatrix}). It is stored in the global option
#' \code{LEAP.verbose} and can also be set directly with
#' \code{options(LEAP.verbose = <level>)}.
#'
#' @param level Verbosity:
#'   \describe{
#'     \item{0}{silent}
#'     \item{1}{main steps with elapsed time (default)}
#'     \item{2}{per-step details (groups, permutations, loaded files)}
#'     \item{3}{debug}
#'   }
#'
#' @returns The level, invisibly.
#' @export
#'
#' @examples
#' fun.SetVerbose(2)   # show per-step details
#' fun.SetVerbose(0)   # silent
#' fun.SetVerbose(1)   # back to the default
fun.SetVerbose <- function(level=1L){
  level <- as.integer(level)
  options(LEAP.verbose=level)
  invisible(level)
}
