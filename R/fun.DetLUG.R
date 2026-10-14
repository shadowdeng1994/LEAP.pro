#' Title
#'
#' @param ObsLUG A dataframe of observing LUGs obtained from fun.LoadInLUGData.
#' @param ShufLUG A dataframe of shuffling LUGs obtained from fun.RunShufLUGEstimator.
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param Threads Number of cores. (retained for compatibility; the computation is vectorised)
#'
#' @returns A dataframe of inferred LUGs.
#' @export
#'
#' @examples Var.LUGs <- fun.DetLUG(Var.Obs,Var.Shuf)
fun.DetLUG <- function(ObsLUG,ShufLUG,p_threshold=0.01,Threads=4){
  ## Empirical P = fraction of permutations (same Total) with Np <= observed Np.
  ## Sorting each Total's Np once lets every gene be scored with a binary
  ## search, so no per-gene filter()/parallel loop is needed.
  tmp.spl <- split(ShufLUG$Np,as.character(ShufLUG$Total))
  tmp.spl <- lapply(tmp.spl,sort)

  tmp.obs <- split(seq_len(nrow(ObsLUG)),as.character(ObsLUG$Total))
  tmp.ep <- rep(NA_real_,nrow(ObsLUG))
  for(tot in names(tmp.obs)){
    s <- tmp.spl[[tot]]
    if(is.null(s)||length(s)==0){ next }
    idx <- tmp.obs[[tot]]
    tmp.ep[idx] <- findInterval(ObsLUG$Np[idx],s)/length(s)
  }

  Var.LUGs <-
    ObsLUG %>%
    mutate(EmpiricalP=tmp.ep) %>%
    mutate(Sig=!is.na(EmpiricalP)&EmpiricalP<p_threshold) %>%
    arrange(EmpiricalP)

  return(Var.LUGs)
}
