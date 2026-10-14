
#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell
#'   phylogeny, or a Seurat object whose cells match the tree tips.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#' @param Assay Assay of the Seurat object to use. (default: NULL, the default assay is used)
#' @param Layer Layer (Seurat >= 5) or slot (Seurat < 5) of the Seurat object to use. (default: "data")
#'
#' @returns NULL.
#' @export
#'
#' @examples fun.RunLUEEstimator(Data.TarCA,Var.TPM)
fun.RunLUGEstimator <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,Threads=4,Redo=FALSE,Assay=NULL,Layer="data"){
  tmp.dirout <- paste0(DirOut,"/LUG/")
  if(!dir.exists(tmp.dirout)){ dir.create(tmp.dirout,recursive = T) }

  tmp.ni <- TarCAObject$NodeIndex
  if(is.null(tmp.ni)){ stop("TarCAObject lacks 'NodeIndex'; recreate it with CreateTarCAObject().") }
  tmp.tipLabels <- tmp.ni$label[tmp.ni$isTip]

  if(inherits(ExpreMatrix,"Seurat")){
    fun.message("--> Extracting expression matrix from Seurat object.")
    ExpreMatrix <- fun.GetExpreMatrix(ExpreMatrix,Assay,Layer)
  }
  if(inherits(ExpreMatrix,"data.frame")){ ExpreMatrix <- as.matrix(ExpreMatrix) }

  ## Align the cells (columns) to the tree tips so the binary annotation can be
  ## used directly against the pre-order index.
  tmp.missing <- setdiff(tmp.tipLabels,colnames(ExpreMatrix))
  if(length(tmp.missing)>0){
    stop("ExpreMatrix is missing ",length(tmp.missing),
         " cell(s) present in TarCAObject$Tree, e.g. ",
         paste(utils::head(tmp.missing,3),collapse=", "),".")
  }
  ExpreMatrix <- ExpreMatrix[,tmp.tipLabels,drop=FALSE]

  if(inherits(ExpreMatrix,"Matrix")){
    tmp.frac <- Matrix::rowSums(ExpreMatrix>0)/ncol(ExpreMatrix)
  }else{
    tmp.frac <- rowMeans(ExpreMatrix>0)
  }
  tmp.expMatrix <- ExpreMatrix[tmp.frac>TopGene,,drop=FALSE]

  tmp.old <- fun.SetPlan(Threads)
  on.exit(future::plan(tmp.old),add=TRUE)

  fun.ParLapply(rownames(tmp.expMatrix),function(ggg){
    tmp.fileout <- paste0(tmp.dirout,"/",ggg,".rds")
    if(!file.exists(tmp.fileout)|Redo){
      tmp.ann <- percent_rank(as.numeric(tmp.expMatrix[ggg,]))>(1-TopGene)
      tmp.res <- fun.LEU_Estimator(tmp.ni,tmp.ann)

      if(is.null(tmp.res)){
        tmp.res <- list(
          PosTips=integer(0),
          FilterBiasParent=data.frame(Filter=logical(0)),
          "LEU_Estimator"=data.frame(TipAnn=TRUE,MonoInfo=NA,Total=sum(tmp.ann),Np=Inf)
        )
      }else{
        ## Only the positive tip positions and the filtered nodes are needed
        ## downstream (avoids storing one label string per cell and gene).
        tmp.res$PosTips <- which(tmp.ann)
      }

      saveRDS(tmp.res,file=tmp.fileout)
      return()
    }
  },Threads)

  message("## Running LUG detector completed.")
}
