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
    fun.message("--> Total: ",sss,"(",match(sss,tmp.set),"/",length(tmp.set),")",level=2)

    tmp.fileout <- paste0(tmp.dirout,"/Total_",sss,".rds")
    tmp.t <- Sys.time()

    if(!file.exists(tmp.fileout)|Redo){
      fun.message("\tRunning...",level=2)

      ## Unique, deterministic seed per (Total, replicate): reproducible
      ## regardless of worker scheduling.
      tmp.totidx <- match(sss,tmp.set)

      tmp.out <-
        fun.ParLapply(1:NumShuf,fun.ShufLUGEstimator_one,Threads,
                      NodeIndex=tmp.ni,ntip=tmp.ntip,sss=sss,
                      SeedFrom=SeedFrom,NumShuf=NumShuf,totidx=tmp.totidx) %>% bind_rows

      saveRDS(tmp.out,file=tmp.fileout)
    }else{
      fun.message("\tLoading: ",tmp.fileout,level=2)
      tmp.out <- readRDS(tmp.fileout)
    }
    fun.message("--> Total ",sss," done in ",fun.secs(tmp.t),"s.",level=2)
    return(tmp.out)
  }) %>% bind_rows

  return(tmp.shuf)
}

#' One permutation replicate (internal worker).
#'
#' @param rrr Replicate index.
#' @param NodeIndex The `NodeIndex` of the tree.
#' @param ntip Number of tips.
#' @param sss Number of TRUE tips for this permutation.
#' @param SeedFrom Base seed (or NULL).
#' @param NumShuf Number of permutations (used to derive the seed).
#' @param totidx Index of the current `Total` (used to derive the seed).
#'
#' @return A one-row data.frame (TipAnn/MonoInfo/Total/Np).
#' @noRd
fun.ShufLUGEstimator_one <- function(rrr,NodeIndex,ntip,sss,SeedFrom,NumShuf,totidx){
  if(!is.null(SeedFrom)){
    ## Pin the RNG kind/sample kind so the stream is identical whether the
    ## workers run sequential or as future multisession processes.
    set.seed(SeedFrom + (totidx-1L)*NumShuf + rrr,
             kind="Mersenne-Twister",sample.kind="Rejection")
  }
  tmp.ann <- sample(1:ntip)<=sss
  tmp.res <- fun.LEU_Estimator(NodeIndex,tmp.ann)
  if(is.null(tmp.res)){
    tmp.res <- data.frame(TipAnn=TRUE,MonoInfo=NA,Total=sum(tmp.ann),Np=Inf)
  }else{
    tmp.res <- tmp.res$LEU_Estimator
  }
  tmp.res
}
