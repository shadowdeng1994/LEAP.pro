#' @title Run label transfer from reference data to query data.
#' @description
#' Run label transfer from reference data to query data, both of which were binarized by LUG.
#'
#' @param QueryMatrix A binarized matrix of query data obtained from GetInternalMatrix.
#' @param RefMatrix A binarized matrix of reference data obtained from GetReferenceMatrix.
#' @param TransferLabel Transfer Labels from reference to query.
#' @param QuerySeurat Existed QuerySeurat.
#' @param RefSeurat Existed RefSeurat.
#' @param Anchors Existed Anchor.
#' @param nPC Number of PC. (default: 30)
#' @param Verbose Verbose or not. (default: FALSE)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param DoUMAP Compute UMAP embedding for the Seurat objects. (default: FALSE)
#'
#' @returns A dataframe of transferred result.
#' @export
#'
#' @examples Var.Predictions <- RunLabelTransfer(Var.InternalMatrix,Var.ReferenceMatrix,Data.Ref_metadata$lineage)
RunLabelTransfer <- function(QueryMatrix=NULL,RefMatrix=NULL,TransferLabel,QuerySeurat=NULL,RefSeurat=NULL,Anchors=NULL,nPC=30,Verbose=FALSE,DirOut="LEAPOutput/",SaveObject=TRUE,DoUMAP=FALSE){
  ## The transfer labels must be in the same cell order as the reference
  ## matrix; if they are named, align them on the reference column names so a
  ## mis-ordered input cannot silently produce wrong labels.
  if(!is.null(RefMatrix)){
    if(is.null(names(TransferLabel))){
      if(length(TransferLabel)!=ncol(RefMatrix)){
        stop("length(TransferLabel) (",length(TransferLabel),") must equal ncol(RefMatrix) (",ncol(RefMatrix),").")
      }
    }else{
      tmp.missing <- setdiff(colnames(RefMatrix),names(TransferLabel))
      if(length(tmp.missing)>0){
        stop("TransferLabel is missing ",length(tmp.missing)," reference cell(s), e.g. ",
             paste(utils::head(tmp.missing,3),collapse=", "),".")
      }
      TransferLabel <- TransferLabel[colnames(RefMatrix)]
    }
  }

  if(SaveObject && !dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }

  if(is.null(Anchors)){
    if(is.null(QuerySeurat)|is.null(RefSeurat)){
      fun.message("--> Creating Seurat object.")
      Var.QuerySeurat <- CreateSeuratObject(counts=QueryMatrix,project = "Query")
      Var.RefSeurat <- CreateSeuratObject(counts=RefMatrix,project = "Ref")
      # ,meta.data = RefMetadata

      fun.message("==> Processing Seurat object.")
      Var.QuerySeurat <- Var.QuerySeurat %>% fun.ProcessSeurat(nPC,Verbose,DoUMAP)
      Var.RefSeurat <- Var.RefSeurat %>% fun.ProcessSeurat(nPC,Verbose,DoUMAP)

      if(SaveObject){
        fun.message("==> Saving Seurat to ",DirOut,".")
        if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
        saveRDS(Var.QuerySeurat,paste0(DirOut,"/Var.QuerySeurat.rds"))
        saveRDS(Var.RefSeurat,paste0(DirOut,"/Var.RefSeurat.rds"))
      }
    }else{
      fun.message("==> Loading Seurat object.")
      Var.QuerySeurat <- QuerySeurat
      Var.RefSeurat <- RefSeurat
    }

    fun.message("==> Finding transfer anchors.")
    tmp.features <- intersect(rownames(Var.QuerySeurat),rownames(Var.RefSeurat))
    if(length(tmp.features)==0){ stop("No shared features between query and reference.") }
    Var.Anchors <- FindTransferAnchors(reference = Var.RefSeurat, query = Var.QuerySeurat, dims = 1:nPC,reference.reduction = "pca",features = tmp.features,verbose=Verbose)

    if(SaveObject){
      fun.message("==> Saving Anchor to ",DirOut,".")
      saveRDS(Var.Anchors,paste0(DirOut,"/Var.Anchors.rds"))
    }
  }else{
    fun.message("==> Loading transfer anchors.")
    Var.Anchors <- Anchors
  }

  fun.message("==> Running label transfer.")
  Var.Predictions <- TransferData(anchorset = Var.Anchors, refdata = TransferLabel, dims = 1:nPC,verbose=Verbose)

  if(SaveObject){
    fun.message("==> Saving Prediction to ",DirOut,".")
    saveRDS(Var.Predictions,paste0(DirOut,"/Var.Predictions.rds"))
  }

  return(Var.Predictions)
}
