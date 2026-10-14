#' Estimate Np for a single binary tip annotation (internal).
#'
#' Equivalent to `TarCA.pro::LEU_Estimator` for one gene/permutation, but
#' rewritten so it scales to large trees:
#'   * subtree counts are obtained with a single vectorised cumulative sum over
#'     the pre-order interval index instead of an O(N) R-level loop;
#'   * the bias-parent filter works on the laminar family of bias parents with
#'     a stack (O(k)) instead of enumerating all O(k^2) pairs.
#'
#' @param NodeIndex The `NodeIndex` of a TarCA ExTree (from CreateTarCAObject).
#' @param tipAnn A logical vector of length `ntip`, in `Tree$tip.label` order.
#'
#' @return NULL if no bias parent is detected, otherwise a list with
#'   `FilterBiasParent` (data.frame with a `Filter` column) and `LEU_Estimator`
#'   (data.frame with TipAnn/MonoInfo/Total/Np).
#' @noRd
fun.LEU_Estimator <- function(NodeIndex,tipAnn){
  ni <- NodeIndex
  N <- ni$N; isTip <- ni$isTip; tin <- ni$tin; tout <- ni$tout; ntip <- ni$ntip
  tipAnn <- as.logical(tipAnn)

  ## Subtree positive counts in one pass: prefix sums over pre-order positions.
  tmp.x <- numeric(N)
  tmp.x[tin[seq_len(ntip)]] <- as.numeric(tipAnn)
  tmp.P <- cumsum(tmp.x)
  tmp.pos <- tmp.P[tout]-c(0,tmp.P)[tin]

  ii <- which(!isTip)
  CladeSize <- ni$size[ii]
  tmp.bg <- sum(tipAnn)/ntip
  Ratio <- tmp.pos[ii]/CladeSize
  padj <- p.adjust(pchisq((tmp.pos[ii]-CladeSize*tmp.bg)^2/(CladeSize*tmp.bg*(1-tmp.bg)),
                          df=1,lower.tail=FALSE),"fdr")
  BiasParent <- (Ratio>tmp.bg) & (padj<0.05) & CladeSize>1
  if(!any(BiasParent)){ return(NULL) }

  bp <- ii[BiasParent]
  tmp.r <- Ratio[BiasParent]
  tinb <- tin[bp]; toutb <- tout[bp]
  k <- length(bp)
  ord <- order(tinb)

  ## nearest bias-parent ancestor (laminar family -> stack), then the maximum
  ## ratio among strict ancestors (`above`) and strict descendants (`below`).
  par <- integer(k); st <- integer(k); top <- 0L
  for(i in ord){
    while(top>0L && toutb[st[top]]<tinb[i]){ top <- top-1L }
    par[i] <- if(top>0L) st[top] else NA_integer_
    top <- top+1L; st[top] <- i
  }
  above <- rep(-Inf,k)
  for(i in ord){ p <- par[i]; if(!is.na(p)){ above[i] <- max(tmp.r[p],above[p]) } }
  below <- rep(-Inf,k)
  for(i in rev(ord)){ p <- par[i]; if(!is.na(p)){ below[p] <- max(below[p],tmp.r[i],below[i]) } }
  Filter <- (tmp.r>above) & (tmp.r>=below)

  tmp.fbp <- data.frame(
    CladeSize=CladeSize[BiasParent],
    Ratio=tmp.r,
    p_chiq0=pchisq((tmp.pos[bp]-CladeSize[BiasParent]*tmp.bg)^2/
                     (CladeSize[BiasParent]*tmp.bg*(1-tmp.bg)),df=1,lower.tail=FALSE),
    padj_chiq0=padj[BiasParent],
    Filter=Filter,
    row.names=ni$label[bp],
    stringsAsFactors=FALSE
  )

  ## Np from the kept (Filter) activation nodes plus the size-1 remainder.
  pf <- tmp.pos[bp][Filter]
  tv <- table(pf)
  cs <- c(1L,as.integer(names(tv)))
  fq <- c(sum(tipAnn)-sum(pf),as.integer(tv))
  o <- order(cs); cs <- cs[o]; fq <- fq[o]
  Total <- as.numeric(sum(cs*fq))

  list(
    FilterBiasParent=tmp.fbp,
    LEU_Estimator=data.frame(
      TipAnn=TRUE,
      MonoInfo=paste0(cs," (",fq,")",collapse=", "),
      Total=Total,
      Np=Total*(Total-1)/sum(cs*(cs-1)*fq)
    )
  )
}
