#' Title
#'
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A dataframe of observed LUGs from local files.
#' @export
#'
#' @examples Var.Obs <- fun.LoadInLUGData()
fun.LoadInLUGData <- function(DirOut="LEAPOutput/",Threads=4){
  tmp.old <- fun.SetPlan(Threads)
  on.exit(future::plan(tmp.old),add=TRUE)

  tmp.ObsLUG <-
    dir(paste0(DirOut,"/LUG/")) %>%
    fun.ParLapply(function(ggg){
      readRDS(paste0(DirOut,"/LUG/",ggg)) %>% .$LEU_Estimator %>% mutate(Symbol=sub("\\.rds$","",ggg))
    },Threads) %>% bind_rows

  return(tmp.ObsLUG)
}
