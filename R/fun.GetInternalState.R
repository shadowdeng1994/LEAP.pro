#' Title
#'
#' @param PhyloDepth A dataframe of PhyloDepth obtained from fun.GetPhyloDepthOfLUGs.
#' @param MinAct Minimum number of activation events for reconstructed internal nodes. (default: 3)
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A matrix of inferred transcriptome of internal nodes.
#' @export
#'
#' @examples fun.GetInternalState(Var.PhyloDepth)
fun.GetInternalState <- function(PhyloDepth,MinAct=3,Threads=4){
  if(is.null(PhyloDepth)||nrow(PhyloDepth)==0){
    return(Matrix::sparseMatrix(i=integer(0),j=integer(0),x=numeric(0),dims=c(0,0)))
  }

  dt <- data.table::as.data.table(PhyloDepth)
  ## One row per (tip, node) and one row per (tip, gene) activation event.
  tmp.nodes <- unique(dt[,list(TipLabel,Node,Node_x=Nor_x)])
  tmp.act   <- unique(dt[,list(TipLabel,Symbol,Symbol_x=Nor_x)])

  ## A gene is active at a node when its activation depth is at or above the
  ## node depth. The non-equi join emits only the TRUE pairs, so the full
  ## per-tip cross product is never materialised (decisive on large trees).
  tmp.true <- tmp.nodes[tmp.act,on=list(TipLabel,Node_x>=Symbol_x),
                        allow.cartesian=TRUE,nomatch=NULL]
  tmp.true <- unique(tmp.true[,list(Node,Symbol)])
  if(nrow(tmp.true)==0){
    return(Matrix::sparseMatrix(i=integer(0),j=integer(0),x=numeric(0),dims=c(0,0)))
  }

  tmp.cnt <- tmp.true[,.N,by=Node]
  tmp.keep <- tmp.cnt[N>=MinAct,Node]
  tmp.true <- tmp.true[Node %in% tmp.keep]

  tmp.genes <- sort(unique(tmp.true$Symbol))
  tmp.nodesOut <- sort(unique(tmp.true$Node))

  Var.InternalState <- Matrix::sparseMatrix(
    i=match(tmp.true$Symbol,tmp.genes),
    j=match(tmp.true$Node,tmp.nodesOut),
    x=1,
    dims=c(length(tmp.genes),length(tmp.nodesOut)),
    dimnames=list(tmp.genes,tmp.nodesOut)
  )

  return(Var.InternalState)
}
