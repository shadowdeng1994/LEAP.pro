#' @title Binarize reference dataset.
#' @description
#' Convert reference dataset into binarized matrix.
#'
#' @param ExpreMatrix A matrix of expression levels corresponding to the reference atlas, or a Seurat object.
#' @param Metadata A dataframe containing metadata for each reference cell. Optional when `ExpreMatrix` is a Seurat object (the object's metadata is used).
#' @param GroupBy Column name used to group reference data. (default: NULL)
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param Threads Number of cores. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param Assay Assay of the Seurat object to use. (default: NULL, the default assay is used)
#' @param Layer Layer (Seurat >= 5) or slot (Seurat < 5) of the Seurat object to use. (default: "data")
#'
#' @returns A binarized sparse matrix (`dgCMatrix`) of reference data.
#' @export
#'
#' @examples GetReferenceMatrix(Data.Ref_Count,Data.Ref_metadata,GroupBy="embryo.time.bin")
#' @examples GetReferenceMatrix(RefSeurat,GroupBy="embryo.time.bin")
GetReferenceMatrix <- function(ExpreMatrix,Metadata=NULL,GroupBy=NULL,TopGene=0.1,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE,Assay=NULL,Layer="data"){
  ## A Seurat object carries both the expression matrix and the cell metadata,
  ## so it can be passed on its own.
  if(inherits(ExpreMatrix,"Seurat")){
    fun.message("--> Extracting expression matrix and metadata from Seurat object.")
    if(is.null(Metadata)){ Metadata <- ExpreMatrix[[]] }
    ExpreMatrix <- fun.GetExpreMatrix(ExpreMatrix,Assay,Layer)
  }
  if(is.null(Metadata)){
    stop("`Metadata` is required when `ExpreMatrix` is not a Seurat object.")
  }

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
    tmp.t <- Sys.time()
    fun.message("--> Group: ",ttt," (",length(tmp.set[[ttt]])," cells)",level=2)

    tmp2 <- ExpreMatrix[,tmp.set[[ttt]],drop=FALSE]
    tmp.mat <- as.matrix(tmp2)
    tmp.n <- ncol(tmp.mat)

    tmp3 <- matrixStats::rowRanks(tmp.mat,ties.method="min")
    tmp.pct <- (tmp3-1)/(tmp.n-1)
    ## a group with a single cell gives 0/0 = NaN; treat it as "not expressed"
    tmp.pct[is.nan(tmp.pct)] <- 0
    ## Keep each binarized group sparse (a gene retains only ~TopGene of the
    ## cells). `Matrix()` on a logical matrix returns a logical sparse matrix;
    ## coercing to "dMatrix" yields a numeric dgCMatrix (1 = retained) that can
    ## be used directly as a Seurat `counts` matrix.
    tmp3 <- methods::as(Matrix::Matrix(tmp.pct > (1-TopGene),sparse=TRUE),"dMatrix")
    dimnames(tmp3) <- dimnames(tmp.mat)

    fun.message("--> Group ",ttt," done in ",fun.secs(tmp.t),"s.",level=2)
    return(tmp3)
  })

  ## every group must keep the same feature order before binding columns
  if(length(tmp.list)>1){
    tmp.rn <- rownames(tmp.list[[1]])
    if(!all(vapply(tmp.list,function(dd) identical(rownames(dd),tmp.rn),logical(1)))){
      stop("Inconsistent feature order across groups.")
    }
  }
  ## `cbind` keeps the group matrices sparse; `dplyr::bind_cols` would coerce
  ## them back to a dense data.frame.
  Var.RefMatrix <- do.call(cbind,tmp.list)

  if(SaveObject){
    fun.message("==> Saving to ",DirOut,".")
    if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
    saveRDS(Var.RefMatrix,paste0(DirOut,"/Var.RefMatrix.rds"))
  }

  return(Var.RefMatrix)
}
