#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param ObsLUG A dataframe of LUGs obtained from fun.LoadInLUGData.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#'
#' @returns A dataframe of shuffling LUGs.
#' @export
#'
#' @examples Var.Shuf <- fun.RunShufLUGEstimator(Data.TarCA,Var.Obs)
fun.RunShufLUGEstimator <- function(TarCAObject,ObsLUG,DirOut="LEAPOutput/",NumShuf=1000,SeedFrom=NULL,Threads=4,Redo=FALSE){
  tmp.dirout <- paste0(DirOut,"/Shuffle/")
  if(!dir.exists(tmp.dirout)){ dir.create(tmp.dirout,recursive = T) }

  tmp.ni <- TarCAObject$NodeIndex
  if(is.null(tmp.ni)){ stop("TarCAObject lacks 'NodeIndex'; recreate it with CreateTarCAObject().") }
  tmp.ntip <- tmp.ni$ntip

  tmp.old <- fun.SetPlan(Threads)
  on.exit(future::plan(tmp.old),add=TRUE)

  tmp.set <- ObsLUG %>% .$Total %>% unique %>% sort
  tmp.shuf <-
  tmp.set %>%
  lapply(function(sss){
    fun.message("--> Total: ",sss,"(",match(sss,tmp.set),"/",length(tmp.set),")")

    tmp.fileout <- paste0(tmp.dirout,"/Total_",sss,".rds")

    if(!file.exists(tmp.fileout)|Redo){
      message("\tRunning...")

      ## Unique, deterministic seed per (Total, replicate): reproducible
      ## regardless of worker scheduling.
      tmp.totidx <- match(sss,tmp.set)

      tmp.out <-
        fun.ParLapply(1:NumShuf,function(rrr){
          if(!is.null(SeedFrom)){
            ## Pin the RNG kind/sample kind so the stream is identical whether
            ## the workers run sequential or as future multisession processes.
            set.seed(SeedFrom + (tmp.totidx-1L)*NumShuf + rrr,
                     kind="Mersenne-Twister",sample.kind="Rejection")
          }
          tmp.ann <- sample(1:tmp.ntip)<=sss
          tmp.res <- fun.LEU_Estimator(tmp.ni,tmp.ann)
          if(is.null(tmp.res)){
            tmp.res <- data.frame(TipAnn=TRUE,MonoInfo=NA,Total=sum(tmp.ann),Np=Inf)
          }else{
            tmp.res <- tmp.res$LEU_Estimator
          }
          return(tmp.res)
        },Threads) %>% bind_rows

      saveRDS(tmp.out,file=tmp.fileout)
    }else{
      message(paste0("\tLoading: ",tmp.fileout))
      tmp.out <- readRDS(tmp.fileout)
    }
    return(tmp.out)
  }) %>% bind_rows

  return(tmp.shuf)
}
