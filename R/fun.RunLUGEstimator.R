
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
#' @param ChunkSize Number of genes per parallel chunk. (default: NULL, i.e.
#'   `ceiling(nGenes / Threads)` -> one chunk per worker)
#'
#' @returns NULL.
#' @export
#'
#' @examples fun.RunLUGEstimator(Data.TarCA,Var.TPM)
fun.RunLUGEstimator <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,Threads=4,Redo=FALSE,Assay=NULL,Layer="data",Force=FALSE,ChunkSize=NULL){
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
    ## Chunk size: default one chunk per worker (ceil(nGenes/Threads)); smaller
    ## chunks with dynamic scheduling (see below) help only when per-gene cost
    ## varies a lot.
    if(is.null(ChunkSize)){ ChunkSize <- ceiling(length(tmp.genes)/max(1L,Threads)) }
    ChunkSize <- max(1L,as.integer(ChunkSize))
    tmp.chunks <- split(tmp.genes,ceiling(seq_along(tmp.genes)/ChunkSize))
    fun.message("--> ",length(tmp.chunks)," chunk(s) of <= ",ChunkSize,
                " genes (Threads=",Threads,").",level=2)

    tmp.old <- fun.SetPlan(Threads)
    on.exit(future::plan(tmp.old),add=TRUE)

    ## Carry each chunk's expression sub-matrix INSIDE X so `future` sends each
    ## worker only its own chunk. Passing the matrix as a global instead makes
    ## future copy the whole matrix to every worker (measured peak ~17x the
    ## data), which is what makes memory blow up.
    tmp.X <- lapply(tmp.chunks,function(gg){
      list(genes=gg,mat=tmp.expMatrix[gg,,drop=FALSE])
    })
    rm(tmp.expMatrix,tmp.chunks)
    gc()

    ## Scheduling=Inf => one future per chunk (dynamic balancing); equal-sized
    ## chunks with one chunk per worker need no balancing, so this is free there.
    fun.ParLapply(tmp.X,fun.LUGEstimator_chunk,Threads,
                  NodeIndex=tmp.ni,TopGene=TopGene,dirout=tmp.dirout,Scheduling=Inf)
    rm(tmp.X)
    gc()
  }

  fun.message("## Running LUG detector completed.")
}

#' Estimate Np for a chunk of genes (internal worker).
#'
#' @param chunk A list with `genes` (character) and `mat` (their expression rows).
#' @param NodeIndex The `NodeIndex` of the tree.
#' @param TopGene Binarisation threshold.
#' @param dirout Directory for the per-gene `.rds`.
#'
#' @return NULL.
#' @noRd
fun.LUGEstimator_chunk <- function(chunk,NodeIndex,TopGene,dirout){
  for(ggg in chunk$genes){
    fun.LUGEstimator_one(ggg,chunk$mat,NodeIndex,TopGene,dirout,TRUE)
  }
  return(NULL)
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
