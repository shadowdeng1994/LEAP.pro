#' @title Calculating potency for each internal node..
#' @description
#' Mapping inferred state to the nodes and their unsolved descendants.
#'
#' The potency itself depends only on the tree and the tip annotation, but it is
#' usually computed for the nodes of one particular state set, so it is always
#' (re)written as `Var.Potency_<index>.rds`, where `<index>` is `Prefix` or a
#' random 4-character one by default.
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param Nodes A vector of internal nodes.
#' @param Threads Number of cores. Kept for compatibility: the computation is
#'   fully vectorised and does not use it. (default: 4)
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param SaveObject Saving all objects. (default: TRUE)
#' @param Prefix Index used to name the output (`Var.Potency_<Prefix>.rds`).
#'   If `NULL`, a random 4-character index of `[0-9A-Za-z]` is generated.
#'   (default: NULL)
#'
#' @returns A dataframe with potency information for each internal node
#'   (`NodeLabel`, `CladeSize`, `On` = occupied categories, `ShannonEntropy` =
#'   plug-in Shannon entropy in nats, `adjSE_MillerMadow` = Miller-Madow
#'   bias-corrected entropy in nats, and `DetailInfo`).
#' @export
#'
#' @examples Var.Potency <- fun.CalculatePotency(Data.TarCA,Var.PredictedState$NodeLabel)
fun.CalculatePotency <- function(TarCAObject,Nodes,Threads=4,DirOut="LEAPOutput/",SaveObject=TRUE,Prefix=NULL){
  tmp.ni <- TarCAObject$NodeIndex
  tmp.tipAnn <- as.character(TarCAObject$Name2Meta$TipAnn)
  names(tmp.tipAnn) <- rownames(TarCAObject$Name2Meta)

  ## Pre-order positions of the tips, grouped by annotation. Counting the
  ## tips of a clade then only needs two binary searches per annotation
  ## (O((N + #nodes*K) log N)) instead of enumerating every clade.
  ##
  ## `findInterval()` requires sorted breakpoints, so sort each annotation's
  ## positions: `NodeIndex$tin` is a DFS visit order, and whether the tips are
  ## visited in tip-label order depends on the child order of the tree's edge
  ## matrix, which does not hold for every tree (on such trees the unsorted
  ## breaks make `findInterval()` fail).
  tmp.tipCat <- tmp.tipAnn[tmp.ni$label[tmp.ni$isTip]]
  tmp.catPos <- lapply(split(tmp.ni$tin[seq_len(tmp.ni$ntip)],tmp.tipCat),sort)
  ## Tips without an annotation cannot be counted; `split()` already drops
  ## them, so drop the `NA` group as well (same as `table()`).
  tmp.cats <- sort(unique(tmp.tipCat))
  tmp.cats <- tmp.cats[!is.na(tmp.cats)]

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

  ## The plug-in entropy above is in bits (log2). The Miller-Madow bias
  ## correction (K-1)/(2N) is defined for the natural-log (nats) estimator, so
  ## convert bits -> nats (x log(2)) first; `adjSE_MillerMadow` is then the
  ## bias-corrected entropy in nats, with K = On (occupied categories) and
  ## N = CladeSize.
  tmp.ent.nats <- tmp.ent*log(2)
  tmp.adjSE <- tmp.ent.nats + (tmp.on-1)/(2*tmp.size)

  tmp.detail <- apply(tmp.counts,1,function(cc){
    paste(paste0(tmp.cats[cc>0],"[",cc[cc>0],"]"),collapse = ";")
  })

  Var.Potency <- data.frame(
    NodeLabel=Nodes,
    CladeSize=tmp.size,
    On=tmp.on,
    ShannonEntropy=tmp.ent.nats,
    adjSE_MillerMadow=tmp.adjSE,
    DetailInfo=tmp.detail,
    stringsAsFactors=FALSE,
    row.names=NULL
  )

  if(SaveObject){
    if(!dir.exists(DirOut)){ dir.create(DirOut,recursive=TRUE) }
    tmp.index <- Prefix
    if(is.null(tmp.index)){ tmp.index <- fun.RandomIndex(4) }
    tmp.fileout <- paste0(DirOut,"/Var.Potency_",tmp.index,".rds")
    fun.message("==> Saving Potency to ",tmp.fileout,".")
    saveRDS(Var.Potency,tmp.fileout)
  }

  return(Var.Potency)
}
