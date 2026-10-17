#' @title Run label transfer from reference data to query data.
#' @description
#' Run label transfer from reference data to query data, both of which were binarized by LUG.
#'
#' The Seurat objects and the transfer anchors do **not** depend on the transfer
#' labels, so they are cached in `DirOut` (`Var.QuerySeurat.rds`,
#' `Var.RefSeurat.rds`, `Var.Anchors.rds`) and reused on the next call (set
#' `Force = TRUE` to recompute them).
#'
#' Only the predictions depend on `TransferLabel`. They are written to
#' `Var.Predictions_<index>.rds`, where `<index>` is `Prefix` or a random
#' 4-character one by default. The prediction file is always (re)computed: a new
#' `Prefix` writes a new file, and re-running with the same `Prefix` overwrites
#' it. Changing only `TransferLabel` therefore lets you run as many label sets as
#' you like without clobbering the others.
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
#' @param Prefix Index used to name the label-dependent output
#'   (`Var.Predictions_<Prefix>.rds`). If `NULL`, a random 4-character index of
#'   `[0-9A-Za-z]` is generated. (default: NULL)
#' @param Force Recompute and overwrite the cached label-independent helpers
#'   `Var.QuerySeurat.rds`, `Var.RefSeurat.rds` and `Var.Anchors.rds` instead of
#'   loading them. (default: FALSE)
#'
#' @returns A dataframe of transferred result.
#' @export
#'
#' @examples Var.Predictions <- RunLabelTransfer(Var.InternalMatrix,Var.ReferenceMatrix,Data.Ref_metadata$lineage)
RunLabelTransfer <- function(QueryMatrix=NULL,RefMatrix=NULL,TransferLabel,QuerySeurat=NULL,RefSeurat=NULL,Anchors=NULL,nPC=30,Verbose=FALSE,DirOut="LEAPOutput/",SaveObject=TRUE,DoUMAP=FALSE,Prefix=NULL,Force=FALSE){
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

  tmp.t0 <- Sys.time()

  ## Fixed names: these only depend on the two matrices, not on the labels.
  tmp.fileQuery <- paste0(DirOut,"/Var.QuerySeurat.rds")
  tmp.fileRef   <- paste0(DirOut,"/Var.RefSeurat.rds")
  tmp.fileAnchor<- paste0(DirOut,"/Var.Anchors.rds")
  ## Label-dependent output: give it an index so different labels can coexist.
  tmp.index <- Prefix
  if(is.null(tmp.index)){ tmp.index <- fun.RandomIndex(4) }
  tmp.filePred <- paste0(DirOut,"/Var.Predictions_",tmp.index,".rds")

  if(is.null(Anchors)){
    if(is.null(QuerySeurat)||is.null(RefSeurat)){
      if(SaveObject&&!Force&&file.exists(tmp.fileQuery)&&file.exists(tmp.fileRef)){
        fun.message("==> Loading cached Seurat objects from ",DirOut,".")
        Var.QuerySeurat <- readRDS(tmp.fileQuery)
        Var.RefSeurat <- readRDS(tmp.fileRef)
      }else{
        fun.message("--> Creating Seurat object.")
        Var.QuerySeurat <- CreateSeuratObject(counts=QueryMatrix,project = "Query")
        Var.RefSeurat <- CreateSeuratObject(counts=RefMatrix,project = "Ref")
        # ,meta.data = RefMetadata

        fun.message("==> Processing Seurat object.")
        Var.QuerySeurat <- Var.QuerySeurat %>% fun.ProcessSeurat(nPC,Verbose,DoUMAP)
        Var.RefSeurat <- Var.RefSeurat %>% fun.ProcessSeurat(nPC,Verbose,DoUMAP)

        if(SaveObject){
          fun.message("==> Saving Seurat to ",DirOut,".")
          saveRDS(Var.QuerySeurat,tmp.fileQuery)
          saveRDS(Var.RefSeurat,tmp.fileRef)
        }
      }
    }else{
      fun.message("==> Loading Seurat object.")
      Var.QuerySeurat <- QuerySeurat
      Var.RefSeurat <- RefSeurat
    }

    if(SaveObject&&!Force&&file.exists(tmp.fileAnchor)){
      fun.message("==> Loading cached anchors from ",DirOut,".")
      Var.Anchors <- readRDS(tmp.fileAnchor)
    }else{
      fun.message("==> Finding transfer anchors.")
      tmp.features <- intersect(rownames(Var.QuerySeurat),rownames(Var.RefSeurat))
      if(length(tmp.features)==0){ stop("No shared features between query and reference.") }
      Var.Anchors <- FindTransferAnchors(reference = Var.RefSeurat, query = Var.QuerySeurat, dims = 1:nPC,reference.reduction = "pca",features = tmp.features,verbose=Verbose)

      if(SaveObject){
        fun.message("==> Saving Anchor to ",DirOut,".")
        saveRDS(Var.Anchors,tmp.fileAnchor)
      }
    }
  }else{
    fun.message("==> Loading transfer anchors.")
    Var.Anchors <- Anchors
  }

  ## The prediction depends on the labels, so it is always (re)computed and
  ## written: a fresh `Prefix` creates a new file, the same `Prefix` overwrites.
  fun.message("==> Running label transfer.")
  Var.Predictions <- TransferData(anchorset = Var.Anchors, refdata = TransferLabel, dims = 1:nPC,verbose=Verbose)

  if(SaveObject){
    fun.message("==> Saving Prediction to ",tmp.filePred,".")
    saveRDS(Var.Predictions,tmp.filePred)
  }

  fun.message("==> RunLabelTransfer total ",fun.secs(tmp.t0),"s.")

  return(Var.Predictions)
}
