#' Title
#'
#' @param ObsLUG A dataframe of observing LUGs obtained from fun.LoadInLUGData.
#' @param ShufLUG A dataframe of shuffling LUGs obtained from fun.RunShufLUGEstimator.
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A dataframe of inferred LUGs.
#' @export
#'
#' @examples Var.LUGs <- fun.DetLUG(Var.Obs,Var.Shuf)
fun.DetLUG <- function(ObsLUG,ShufLUG,p_threshold=0.01,Threads=4){
  Var.LUGs <-
  1:nrow(ObsLUG) %>%
  mclapply(mc.cores = Threads,function(iii){
    tmp1 <- ObsLUG[iii,]
    tmp2 <- ShufLUG %>% filter(Total==ObsLUG$Total[iii]) %>% .$Np
    tmp1 %>% mutate(EmpiricalP=if(length(tmp2)==0){ NA_real_ }else{ mean(tmp2<=Np) })
  }) %>% bind_rows %>%
  mutate(Sig=!is.na(EmpiricalP)&EmpiricalP<p_threshold) %>%
  arrange(EmpiricalP)

  return(Var.LUGs)
}
