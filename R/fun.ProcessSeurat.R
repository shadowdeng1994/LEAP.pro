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

  ## On its parallel branch (`future::nbrOfWorkers() > 1`) `NormalizeData`
  ## calls `Seurat:::Parenting(parent.find = "Seurat")`, which warns unless
  ## exactly one stack frame's call text contains "Seurat". Our frames
  ## ("fun.ProcessSeurat" and the pipe text "SeuratObject %>% NormalizeData")
  ## add two more matches, so Seurat gives up and warns. Nothing is lost: the
  ## values it would inject are local to that call, and `data`, `scale.data` and
  ## the PCA embedding are identical with and without the warning. Muffle
  ## exactly that message and keep every other warning visible.
  tmp.muffle <- function(w){
    if(identical(conditionMessage(w),"Cannot find a parent environment called Seurat")){
      invokeRestart("muffleWarning")
    }
  }

  SeuratObject <-
    withCallingHandlers({
      SeuratObject %>%
      NormalizeData(verbose = Verbose) %>%
      ScaleData(verbose = Verbose) %>%
      RunPCA(features = rownames(.),npcs=tmp.npcs,verbose = Verbose)
    },warning = tmp.muffle)

  ## UMAP is only used for visualisation and is not required by the label
  ## transfer, so it is skipped by default to save substantial time. Clamp the
  ## dimensions to the number of PCs actually computed (few features -> error).
  if(DoUMAP){
    SeuratObject <- RunUMAP(SeuratObject,dims=1:min(nPC,tmp.npcs),verbose=Verbose)
  }

  return(SeuratObject)
}
