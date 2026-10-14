#' Title
#'
#' @param TarCAObject A TarCA object obtained from CreateTarCAObject.
#' @param LUG A dataframe of LUGs obtained from fun.ExtractLUG.
#' @param DirOut Path to saving the LUG object. (default: "LEAPOutput/")
#' @param Threads Number of cores. (default: 4)
#'
#' @returns A dataframe of normalized depth of internal nodes with LUG activation events.
#' @export
#'
#' @examples fun.GetPhyloDepthOfLUGs(Data.TarCA,Var.LUGs)
fun.GetPhyloDepthOfLUGs <- function(TarCAObject,LUG,DirOut="LEAPOutput/",Threads=4){
  tmp.ni <- TarCAObject$NodeIndex
  ## Nor_x is the same for every gene, so compute it once.
  tmp.Nor_x <- fun.GetNodeDepth(TarCAObject$Tree)
  tmp.Nor_x <- tmp.Nor_x/median(tmp.Nor_x[tmp.ni$label[tmp.ni$isTip]])

  ## Named lookups avoid per-gene joins.
  tmp.cellAnn <- TarCAObject$Ann$TipAnn
  names(tmp.cellAnn) <- TarCAObject$Ann$TipLabel
  tmp.tin <- tmp.ni$tin; names(tmp.tin) <- tmp.ni$label
  tmp.tout <- tmp.ni$tout; names(tmp.tout) <- tmp.ni$label
  tmp.tipLabels <- tmp.ni$label[tmp.ni$isTip]

  tmp.old <- fun.SetPlan(Threads)
  on.exit(future::plan(tmp.old),add=TRUE)

  Var.PhyloDepthOfLUGs <-
  LUG %>% filter(Sig) %>% .$Symbol %>%
  fun.ParLapply(fun.PhyloDepth_one,Threads,
                DirOut=DirOut,NodeIndex=tmp.ni,Nor_x=tmp.Nor_x,cellAnn=tmp.cellAnn,
                tin=tmp.tin,tout=tmp.tout,tipLabels=tmp.tipLabels)

  Var.PhyloDepthOfLUGs <- as.data.frame(data.table::rbindlist(Var.PhyloDepthOfLUGs,use.names=TRUE,fill=TRUE))

  return(Var.PhyloDepthOfLUGs)
}

#' Map one gene's activation nodes to the depths of their positive tips (internal worker).
#'
#' @param ggg Gene symbol.
#' @param DirOut Directory holding `LUG/<gene>.rds`.
#' @param NodeIndex The `NodeIndex` of the tree.
#' @param Nor_x Normalised node depths, named by node label.
#' @param cellAnn Cell-type annotation, named by tip label.
#' @param tin,tout Named pre-order in/out positions.
#' @param tipLabels Tip labels in tree order.
#'
#' @return A data.frame (TipLabel/TipAnn/Node/Nor_x/Symbol) or NULL.
#' @noRd
fun.PhyloDepth_one <- function(ggg,DirOut,NodeIndex,Nor_x,cellAnn,tin,tout,tipLabels){
  tmp1 <- readRDS(paste0(DirOut,"/LUG/",ggg,".rds"))

  tmp.fbp <- tmp1$FilterBiasParent
  if(is.null(tmp.fbp)||nrow(tmp.fbp)==0){ return(NULL) }
  tmp.nodes <- rownames(tmp.fbp)[tmp.fbp$Filter]
  if(length(tmp.nodes)==0){ return(NULL) }

  if(length(tmp1$PosTips)==0){ return(NULL) }
  tmp.posTips <- tipLabels[tmp1$PosTips]

  ## Index the positive tips by pre-order position once, then slice the band
  ## that falls inside each filtered node with two binary searches. Filtered
  ## nodes may nest, so a tip can belong to several nodes (as in the original
  ## implementation); cost is O(#nodes*log n + output).
  tmp.pos <- tin[tmp.posTips]
  tmp.o <- order(tmp.pos)
  tmp.posS <- tmp.pos[tmp.o]
  tmp.tipS <- tmp.posTips[tmp.o]

  tmp.res <- vector("list",length(tmp.nodes))
  for(i in seq_along(tmp.nodes)){
    nnn <- tmp.nodes[i]
    lo <- findInterval(tin[[nnn]]-1L,tmp.posS)+1L
    hi <- findInterval(tout[[nnn]],tmp.posS)
    if(hi>=lo){
      tmp.res[[i]] <- data.frame(TipLabel=tmp.tipS[lo:hi],Node=nnn,stringsAsFactors=FALSE)
    }
  }
  tmp3 <- data.table::rbindlist(tmp.res)
  if(nrow(tmp3)==0){ return(NULL) }

  data.frame(
    TipLabel=tmp3$TipLabel,
    TipAnn=as.character(cellAnn[tmp3$TipLabel]),
    Node=tmp3$Node,
    Nor_x=as.numeric(Nor_x[tmp3$Node]),
    Symbol=ggg,
    stringsAsFactors=FALSE
  )
}
