
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
#' @param Force Recompute genes whose `LUG/*.rds` already exists. (default: FALSE)
#'
#' @returns NULL.
#' @export
#'
#' @examples fun.RunLUEEstimator(Data.TarCA,Var.TPM)
fun.RunLUGEstimator <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,Threads=4,Redo=FALSE,Assay=NULL,Layer="data",Force=FALSE){
  Redo <- isTRUE(Redo)||isTRUE(Force)
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

  ## Release the full expression matrix and the filter vector early: otherwise
  ## they stay alive in this frame for the whole (long) gene loop.
  rm(ExpreMatrix,tmp.frac)
  gc()

  ## Only (re)compute the genes that still need it; this also avoids exporting
  ## the expression matrix when there is nothing to do.
  tmp.genes <- rownames(tmp.expMatrix)
  if(!Redo){ tmp.genes <- tmp.genes[!file.exists(paste0(tmp.dirout,"/",tmp.genes,".rds"))] }
  fun.message("--> ",length(tmp.genes)," genes to compute (",
              nrow(tmp.expMatrix)-length(tmp.genes)," cached).",level=2)

  if(length(tmp.genes)>0){
    tmp.todo <- tmp.expMatrix[tmp.genes,,drop=FALSE]
    rm(tmp.expMatrix)
    gc()

    tmp.old <- fun.SetPlan(Threads)
    on.exit(future::plan(tmp.old),add=TRUE)

    ## The worker is a package-level function so `future` serialises only it and
    ## the objects passed explicitly, instead of the whole calling frame (which
    ## holds the full expression matrix / Seurat / TarCA object).
    fun.ParLapply(tmp.genes,fun.LUGEstimator_one,Threads,
                  expMatrix=tmp.todo,NodeIndex=tmp.ni,TopGene=TopGene,
                  dirout=tmp.dirout,Redo=TRUE)
    rm(tmp.todo)
    gc()
  }

  fun.message("## Running LUG detector completed.")
}

#' Estimate Np for a single gene (internal worker).
#'
#' @param ggg Gene name (row of `expMatrix`).
#' @param expMatrix Genes x cells expression matrix (tree-tip columns).
#' @param NodeIndex The `NodeIndex` of the tree.
#' @param TopGene Binarisation threshold.
#' @param dirout Directory for the per-gene `.rds`.
#' @param Redo Recompute even if the file exists.
#'
#' @return NULL (writes a file).
#' @noRd
fun.LUGEstimator_one <- function(ggg,expMatrix,NodeIndex,TopGene,dirout,Redo){
  tmp.fileout <- paste0(dirout,"/",ggg,".rds")
  if(!file.exists(tmp.fileout)|Redo){
    tmp.ann <- percent_rank(as.numeric(expMatrix[ggg,]))>(1-TopGene)
    tmp.res <- fun.LEU_Estimator(NodeIndex,tmp.ann)

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
  }
  return(NULL)
}
