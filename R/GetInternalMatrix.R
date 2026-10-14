#' @title Reconstruct the LUG-transcriptome of internal nodes.
#' @description
#' Reconstruct the transcriptome of internal nodes based on LUG inference.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject. (default: NULL)
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell phylogeny, or a Seurat object whose cells match the tree tips. (default: NULL)
#' @param LUGs Either a dataframe of LUGs obtained from fun.ExtractLUG, or a plain character vector of LUG gene symbols. If supplied, the LUG-estimation step is skipped (the per-gene `LUG/*.rds` files must already exist in `DirOut`). (optional,default: NULL)
#' @param PhyloDepth A dataframe of PhyloDepth obtained from fun.GetPhyloDepthOfLUGs. (optional,default: NULL)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param MinAct Minimum number of activation events for reconstructed internal nodes. (default: 3)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param Redo Redo LUGs searching. (default: FALSE)
#' @param Assay Assay of the Seurat object to use. (default: NULL, the default assay is used)
#' @param Layer Layer (Seurat >= 5) or slot (Seurat < 5) of the Seurat object to use. (default: "data")
#' @param OnlyLUG If TRUE, run only the LUG-extraction step and return the LUGs (skip PhyloDepth conversion and internal-state inference). (default: FALSE)
#'
#' @returns A matrix of inferred transcriptome of internal nodes, or a dataframe of LUGs when `OnlyLUG = TRUE`.
#' @export
#'
#' @examples Data.TarCA <- CreateTarCAObject(Var.Tree,Var.Ann)
#' @examples Var.InternalState <- GetInternalMatrix(Data.TarCA,Var.TPM)
#' @examples Var.LUGs <- GetInternalMatrix(Data.TarCA,Var.TPM,OnlyLUG=TRUE)
#' @examples Var.InternalState <- GetInternalMatrix(Data.TarCA,LUGs=c("GeneA","GeneB"))
GetInternalMatrix <- function(TarCAObject=NULL,ExpreMatrix=NULL,LUGs=NULL,PhyloDepth=NULL,DirOut="LEAPOutput/",TopGene=0.1,NumShuf=1000,p_threshold=0.01,MinAct=3,SeedFrom=NULL,Threads=4,SaveObject=TRUE,Redo=FALSE,Assay=NULL,Layer="data",OnlyLUG=FALSE){
  tmp.t0 <- Sys.time()

  ## accept a plain character vector of gene symbols as the LUG set
  if(!is.null(LUGs)&&!is.data.frame(LUGs)){
    Var.LUGs <- data.frame(Symbol=as.character(LUGs),Sig=TRUE,stringsAsFactors=FALSE)
  }else{
    Var.LUGs <- LUGs
  }

  if(is.null(Var.LUGs)){
    fun.message("==> Extracting LUG.")
    tmp.t <- Sys.time()
    Var.LUGs <- fun.ExtractLUG(TarCAObject,ExpreMatrix,DirOut,TopGene,NumShuf,p_threshold,SeedFrom,Threads,Redo,Assay,Layer)
    fun.message("==> Extracting LUG done in ",fun.secs(tmp.t),"s.")
  }else{
    fun.message("==> Using supplied LUGs (",nrow(Var.LUGs)," genes); skipping LUG estimation.")
  }

  ## only run the LUG step and return it
  if(OnlyLUG){
    if(SaveObject){
      if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
      saveRDS(Var.LUGs,paste0(DirOut,"/Var.LUGs.rds"))
    }
    fun.message("==> OnlyLUG = TRUE; returning LUGs (total ",fun.secs(tmp.t0),"s).")
    return(Var.LUGs)
  }

  if(is.null(PhyloDepth)){
    fun.message("==> Converting to PhyloDepth.")
    tmp.t <- Sys.time()
    Var.PhyloDepthOfLUGs <- fun.GetPhyloDepthOfLUGs(TarCAObject,Var.LUGs,DirOut=DirOut,Threads)
    fun.message("==> Converting to PhyloDepth done in ",fun.secs(tmp.t),"s.")
  }else{
    fun.message("==> Using supplied PhyloDepth.")
    Var.PhyloDepthOfLUGs <- PhyloDepth
  }

  fun.message("==> Inferring internal state.")
  tmp.t <- Sys.time()
  Var.InternalState <- fun.GetInternalState(Var.PhyloDepthOfLUGs,MinAct,Threads)
  fun.message("==> Inferring internal state done in ",fun.secs(tmp.t),"s.")

  if(SaveObject){
    fun.message("==> Saving to ",DirOut,".")
    tmp.t <- Sys.time()
    saveRDS(Var.LUGs,paste0(DirOut,"/Var.LUGs.rds"))
    saveRDS(Var.PhyloDepthOfLUGs,paste0(DirOut,"/Var.PhyloDepthOfLUGs.rds"))
    saveRDS(Var.InternalState,paste0(DirOut,"/Var.InternalState.rds"))
    fun.message("==> Saving done in ",fun.secs(tmp.t),"s.")
  }

  fun.message("==> GetInternalMatrix total ",fun.secs(tmp.t0),"s.")

  return(Var.InternalState)
}
