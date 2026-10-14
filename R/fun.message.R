#' Title
#'
#' @param ... Message you want to print.
#' @param level Minimum verbosity level (`options(LEAP.verbose)`, default 1) at
#'   which this message is printed. See \code{\link{fun.SetVerbose}}.
#'
#' @returns Message with time.
#' @export
#'
#' @examples fun.message("Hello.")
fun.message <- function(...,level=1){
  if(getOption("LEAP.verbose",1L)>=level){
    cat(paste0("[",date(),"]"),...)
    cat("\n")
    flush.console()
  }
}

#' Elapsed seconds since a start time (internal).
#'
#' @param t0 A `POSIXct` start time (e.g. `Sys.time()`).
#'
#' @return Numeric seconds, rounded to 2 decimals.
#' @noRd
fun.secs <- function(t0){
  round(as.numeric(difftime(Sys.time(),t0,units="secs")),2)
}
