#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Predictions A dataframe of transferred result.
#'
#' @returns A dataframe with later state.
#' @export
#'
#' @examples fun.AssignLaterState(Data.TarCA,Var.Predictions)
fun.AssignLaterState <- function(TarCAObject,Predictions){
  tmp.ni <- TarCAObject$NodeIndex
  if(is.null(tmp.ni)){ stop("TarCAObject lacks 'NodeIndex'; recreate it with CreateTarCAObject().") }

  tmp.label <- tmp.ni$label
  tmp.parent <- tmp.ni$parent
  tmp.isTip <- tmp.ni$isTip

  ## Iterative pre-order (parent before child) replaces the recursive
  ## traversal: no C-stack limit on deep trees and no per-node rbind (O(N^2)).
  tmp.order <- tmp.ni$preorder[!tmp.isTip[tmp.ni$preorder]]

  ## Exact row-name matching. The original used `Predictions[NNNode, ]`, which
  ## performs data.frame partial matching (e.g. "Node_20" wrongly matching
  ## "Node_200"); exact matching fixes that.
  tmp.j <- match(tmp.label,rownames(Predictions))
  tmp.pred <- setNames(rep(NA_character_,length(tmp.label)),tmp.label)
  tmp.ok <- !is.na(tmp.j)
  tmp.pred[tmp.ok] <- as.character(Predictions[["TransLabel"]])[tmp.j[tmp.ok]]

  tmp.state <- character(length(tmp.label))
  tmp.group <- character(length(tmp.label))
  tmp.anc <- character(length(tmp.label))

  for(v in tmp.order){
    lab <- tmp.label[v]
    tmp.cur <- tmp.pred[[lab]]
    if(is.na(tmp.cur)){
      if(tmp.parent[v]==v){          ## root with no prediction
        tmp.state[v] <- "Root"
        tmp.group[v] <- "Legacy"
        tmp.anc[v] <- lab
      }else{
        tmp.state[v] <- tmp.state[tmp.parent[v]]
        tmp.group[v] <- "Legacy"
        tmp.anc[v] <- tmp.anc[tmp.parent[v]]
      }
    }else{
      tmp.state[v] <- tmp.cur
      tmp.group[v] <- "LUG"
      tmp.anc[v] <- lab
    }
  }

  data.frame(
    NodeLabel=tmp.label[tmp.order],
    TransLabel=tmp.state[tmp.order],
    Type=tmp.group[tmp.order],
    AncLabel=tmp.anc[tmp.order],
    stringsAsFactors=FALSE
  )
}
