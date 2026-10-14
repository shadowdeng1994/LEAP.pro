#' @title Binarize reference dataset.
#' @description
#' Convert reference dataset into binarized matrix.
#'
#' @param ExpreMatrix A matrix of expression levels corresponding to the reference atlas.
#' @param Metadata A dataframe containing metadata for each reference cell.
#' @param GroupBy Column name used to group reference data. (default: NULL)
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param Threads Number of cores. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A binarized matrix of reference data.
#' @export
#'
#' @examples GetReferenceMatrix(Data.Ref_Count,Data.Ref_metadata,GroupBy="embryo.time.bin")
GetReferenceMatrix <- function(ExpreMatrix,Metadata,GroupBy=NULL,TopGene=0.1,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE){
  if(!is.null(GroupBy)){
    tmp.set <- split(x=rownames(Metadata),f=Metadata[,GroupBy])
  }else{
    tmp.set <- list("All"=rownames(Metadata))
  }

  ## Ranking is done group-wise with a single vectorised C++ call per group
  ## (`matrixStats::rowRanks`, ties.method="min" == dplyr::percent_rank)
  ## instead of one R-level percent_rank per gene.
  tmp.list <-
  names(tmp.set) %>%
  lapply(function(ttt){
    fun.message("--> Group: ",ttt)

    tmp2 <- ExpreMatrix[,tmp.set[[ttt]],drop=FALSE]
    tmp.mat <- as.matrix(tmp2)
    tmp.n <- ncol(tmp.mat)

    tmp3 <- matrixStats::rowRanks(tmp.mat,ties.method="min")
    tmp.pct <- (tmp3-1)/(tmp.n-1)
    ## a group with a single cell gives 0/0 = NaN; treat it as "not expressed"
    tmp.pct[is.nan(tmp.pct)] <- 0
    tmp3 <- as.data.frame((tmp.pct > (1-TopGene))*1)
    dimnames(tmp3) <- dimnames(tmp.mat)

    return(tmp3)
  })

  ## every group must keep the same feature order before binding columns
  if(length(tmp.list)>1){
    tmp.rn <- rownames(tmp.list[[1]])
    if(!all(vapply(tmp.list,function(dd) identical(rownames(dd),tmp.rn),logical(1)))){
      stop("Inconsistent feature order across groups.")
    }
  }
  Var.RefMatrix <- bind_cols(tmp.list)

  if(SaveObject){
    fun.message("==> Saving to ",DirOut,".")
    if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
    saveRDS(Var.RefMatrix,paste0(DirOut,"/Var.RefMatrix.rds"))
  }

  return(Var.RefMatrix)
}
