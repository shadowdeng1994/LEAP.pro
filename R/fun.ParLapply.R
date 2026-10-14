#' Choose a parallel backend for a given number of workers (internal).
#'
#' Returns the current `future::plan()` and, when `Threads > 1` and the current
#' plan is sequential, switches to a multisession backend with `Threads`
#' workers. The caller is expected to restore the previous plan with
#' `on.exit(future::plan(old), add = TRUE)`. A non-sequential plan set by the
#' user is respected.
#'
#' @param Threads Number of workers. (default: 1)
#'
#' @return The previous `future` plan.
#' @noRd
fun.SetPlan <- function(Threads=1){
  old <- future::plan()
  if(!is.null(Threads) && Threads>1 && inherits(old,"sequential")){
    future::plan(future::multisession,workers=Threads)
  }
  old
}

#' Parallel lapply honouring `Threads` (internal).
#'
#' Uses `future.apply::future_lapply` (multisession), which works on every
#' platform (including Windows) and avoids the fork-related failures of
#' `parallel::mclapply` inside RStudio. The plan should already be configured
#' by `fun.SetPlan`.
#'
#' @param X A vector or list.
#' @param FUN Function to apply to each element.
#' @param Threads Number of workers. (default: 1)
#' @param ... Further arguments passed to `FUN`.
#' @param Scheduling Passed to `future.apply::future_lapply` as
#'   `future.scheduling`; `Inf` schedules one future per element (dynamic load
#'   balancing, useful when element costs vary). (default: NULL, future's own
#'   default)
#'
#' @return A list.
#' @noRd
fun.ParLapply <- function(X,FUN,Threads=1,...,Scheduling=NULL){
  if(is.null(Threads)||Threads<=1||length(X)<=1){ return(lapply(X,FUN,...)) }
  if(is.null(Scheduling)){
    future.apply::future_lapply(X,FUN,...,future.seed=TRUE)
  }else{
    future.apply::future_lapply(X,FUN,...,future.seed=TRUE,future.scheduling=Scheduling)
  }
}
