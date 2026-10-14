#' Extract a genes x cells expression matrix (internal).
#'
#' Accepts either a Seurat object or an ordinary matrix/data.frame and always
#' returns a genes x cells matrix (genes as rows, cells as columns) ready for
#' the LUG estimator. For a Seurat object the requested assay/layer is used.
#'
#' @param x A Seurat object, matrix or data.frame.
#' @param Assay Assay to extract from a Seurat object. (default: NULL, the
#'   default assay is used)
#' @param Layer Layer (Seurat >= 5) or slot (Seurat < 5) to extract from a
#'   Seurat object. (default: "data")
#'
#' @return A genes x cells expression matrix.
#' @noRd
fun.GetExpreMatrix <- function(x,Assay=NULL,Layer="data"){
  if(!inherits(x,"Seurat")){ return(x) }
  if(is.null(Assay)){ Assay <- Seurat::DefaultAssay(x) }
  tryCatch(
    Seurat::GetAssayData(x,assay=Assay,layer=Layer),
    error=function(e) Seurat::GetAssayData(x,assay=Assay,slot=Layer)
  )
}
