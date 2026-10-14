#' Title
#'
#' @param SeuratObject A Seurat object.
#' @param nPC Number of PC.
#' @param Verbose Verbose or not.
#' @param DoUMAP Compute UMAP embedding. (default: FALSE)
#'
#' @returns A processed Seurat object
#' @export
#'
#' @examples fun.ProcessSeurat(Var.QuerySeurat)
fun.ProcessSeurat <- function(SeuratObject,nPC=30,Verbose=FALSE,DoUMAP=FALSE){
  tmp.npcs <- min(nrow(SeuratObject),100)

  SeuratObject <-
    SeuratObject %>%
    NormalizeData(verbose = Verbose) %>%
    ScaleData(verbose = Verbose) %>%
    RunPCA(features = rownames(.),npcs=tmp.npcs,verbose = Verbose)

  ## UMAP is only used for visualisation and is not required by the label
  ## transfer, so it is skipped by default to save substantial time. Clamp the
  ## dimensions to the number of PCs actually computed (few features -> error).
  if(DoUMAP){
    SeuratObject <- RunUMAP(SeuratObject,dims=1:min(nPC,tmp.npcs),verbose=Verbose)
  }

  return(SeuratObject)
}
