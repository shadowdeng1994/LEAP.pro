#' @title Extract the LUGs from cell phylogeny and corresponding transcritptome.
#' @description
#' LUGs inference based on permutation test.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param ExpreMatrix A matrix of expression levels corresponding to the cell phylogeny.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param TopGene Threshold of highly expressed population. (default: 0.1)
#' @param NumShuf Number of shuffling. (default: 1000)
#' @param p_threshold Threshold of empirical P-value. (default: 0.01)
#' @param SeedFrom Set the beginning of seed. (default: NULL)
#' @param Threads Number of cores. (default: 4)
#' @param Redo Redo LUGs searching. (default: FALSE)
#' @param Assay Assay of the Seurat object to use. (default: NULL, the default assay is used)
#' @param Layer Layer (Seurat >= 5) or slot (Seurat < 5) of the Seurat object to use. (default: "data")
#'
#' @returns A dataframe of inferred LUGs.
#' @export
#'
#' @examples Data.TarCA <- CreateTarCAObject(Var.Tree,Var.Ann)
#' @examples Var.LUGs <- fun.ExtractLUG(Data.TarCA,Var.TPM)
fun.ExtractLUG <- function(TarCAObject,ExpreMatrix,DirOut="LEAPOutput/",TopGene=0.1,NumShuf=1000,p_threshold=0.01,SeedFrom=NULL,Threads=4,Redo=FALSE,Assay=NULL,Layer="data"){
  if(!dir.exists(DirOut)){
    fun.message("--> ",DirOut,"is absent. Creating.")
    dir.create(DirOut)
  }

  fun.message("--> Running LUG estimator.")
  tmp.t <- Sys.time()
  fun.RunLUGEstimator(TarCAObject,ExpreMatrix,DirOut,TopGene,Threads,Redo,Assay,Layer)
  fun.message("--> LUG estimator done in ",fun.secs(tmp.t),"s.")

  tmp.t <- Sys.time()
  Var.Obs <- fun.LoadInLUGData(DirOut,Threads)
  fun.message("--> Loaded ",nrow(Var.Obs)," observed LUGs in ",fun.secs(tmp.t),"s.",level=2)

  fun.message("--> Shuffling (",NumShuf,"times).")
  tmp.t <- Sys.time()
  Var.Shuf <- fun.RunShufLUGEstimator(TarCAObject,Var.Obs,DirOut,NumShuf,SeedFrom,Threads,Redo)
  fun.message("--> Shuffling done in ",fun.secs(tmp.t),"s.")

  fun.message("--> Extracting high confident LUGs.")
  tmp.t <- Sys.time()
  Var.LUGs <- fun.DetLUG(Var.Obs,Var.Shuf,p_threshold,Threads)
  fun.message("--> LUG detection done in ",fun.secs(tmp.t),"s. (",sum(Var.LUGs$Sig)," significant)")
  return(Var.LUGs)
}
