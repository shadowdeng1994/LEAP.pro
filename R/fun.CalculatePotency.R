#' @title Calculating potency for each internal node..
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Nodes A vector of internal nodes.
#' @param Threads Number of cores. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#'
#' @returns A dataframe with potency information for each internal node.
#' @export
#'
#' @examples Var.Potency <- CalculatePotency(Data.TarCA,Var.PredictedState$NodeLabel)
fun.CalculatePotency <- function(TarCAObject,Nodes,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE){
  tmp.ni <- TarCAObject$NodeIndex
  tmp.tipAnn <- as.character(TarCAObject$Name2Meta$TipAnn)
  names(tmp.tipAnn) <- rownames(TarCAObject$Name2Meta)

  ## Pre-order positions of the tips, grouped by annotation. Counting the
  ## tips of a clade then only needs two binary searches per annotation
  ## (O((N + #nodes*K) log N)) instead of enumerating every clade.
  tmp.tipCat <- tmp.tipAnn[tmp.ni$label[tmp.ni$isTip]]
  tmp.catPos <- split(tmp.ni$tin[seq_len(tmp.ni$ntip)],tmp.tipCat)
  tmp.cats <- sort(unique(tmp.tipCat))

  tmp.idx <- match(Nodes,tmp.ni$label)
  tmp.tins <- tmp.ni$tin[tmp.idx]
  tmp.touts <- tmp.ni$tout[tmp.idx]

  tmp.counts <- vapply(tmp.cats,function(cc){
    tmp.pa <- tmp.catPos[[cc]]
    findInterval(tmp.touts,tmp.pa)-findInterval(tmp.tins-1L,tmp.pa)
  },numeric(length(Nodes)))
  if(is.null(dim(tmp.counts))){ tmp.counts <- matrix(tmp.counts,nrow=length(Nodes)) }
  colnames(tmp.counts) <- tmp.cats

  tmp.size <- rowSums(tmp.counts)
  tmp.on <- rowSums(tmp.counts>0)
  tmp.denom <- tmp.size; tmp.denom[tmp.denom==0] <- 1
  tmp.p <- tmp.counts/tmp.denom
  tmp.logp <- tmp.p*log2(tmp.p)
  tmp.logp[tmp.counts==0] <- 0
  tmp.ent <- -rowSums(tmp.logp)

  tmp.detail <- apply(tmp.counts,1,function(cc){
    paste(paste0(tmp.cats[cc>0],"[",cc[cc>0],"]"),collapse = ";")
  })

  Var.Potency <- data.frame(
    NodeName=Nodes,
    CladeSize=tmp.size,
    On=tmp.on,
    ShannonEntropy=tmp.ent,
    DetailInfo=tmp.detail,
    stringsAsFactors=FALSE,
    row.names=NULL
  )

  if(SaveObject){
    fun.message("==> Saving Potency to ",DirOut,".")
    saveRDS(Var.Potency,paste0(DirOut,"/Var.Potency.rds"))
  }

  return(Var.Potency)
}
