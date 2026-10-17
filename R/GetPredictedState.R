#' @title Get Inferred states from label transfer data.
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' The result depends on `Predictions` (i.e. on the transferred labels), so it
#' is always (re)computed and written to `Var.PredictedState_<index>.rds`, where
#' `<index>` is `Prefix` or a random 4-character one by default. A new `Prefix`
#' writes a new file, and re-running with the same `Prefix` overwrites it.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Predictions A dataframe of transferred result.
#' @param RemoveRoot Remove root's name.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param Prefix Index used to name the output (`Var.PredictedState_<Prefix>.rds`).
#'   If `NULL`, a random 4-character index of `[0-9A-Za-z]` is generated.
#'   (default: NULL)
#'
#' @returns A datafraome of impuated annotation on the cell phylogeny.
#' @export
#'
#' @examples Var.PredictedState <- GetPredictedState(Data.TarCA,Var.Predictions)
GetPredictedState <- function(TarCAObject,Predictions,RemoveRoot=NULL,DirOut="LEAPOutput/",SaveObject=TRUE,Prefix=NULL){
  Var.PredictedState <-
  fun.AssignLaterState(
    TarCAObject,
    Predictions %>% select(TransLabel=predicted.id)
  )

  ## Map each node's nearest ancestral LUG node (AncLabel) to its predicted
  ## label with exact matching. The previous `Predictions[AncLabel, ]` used
  ## data.frame partial matching and could silently drop rows when it
  ## returned NA.
  tmp.idx <- match(Var.PredictedState$AncLabel,rownames(Predictions))
  Var.PredictedState$TransLabel <- Predictions[["predicted.id"]][tmp.idx]

  ## Drop nodes that cannot be assigned to any ancestral state (no predicted
  ## ancestor, i.e. `NA`) and, optionally, a specific label such as "Root".
  if(!is.null(RemoveRoot)){
    Var.PredictedState <- Var.PredictedState %>%
      filter(!is.na(TransLabel),TransLabel!=RemoveRoot)
  }

  if(SaveObject){
    if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
    tmp.index <- Prefix
    if(is.null(tmp.index)){ tmp.index <- fun.RandomIndex(4) }
    tmp.fileout <- paste0(DirOut,"/Var.PredictedState_",tmp.index,".rds")
    fun.message("==> Saving PredictedState to ",tmp.fileout,".")
    saveRDS(Var.PredictedState,tmp.fileout)
  }

  return(Var.PredictedState)
}
