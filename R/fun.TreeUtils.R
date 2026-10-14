#' Compute root-to-node distances named by node label (internal).
#'
#' Reproduces the `x` column that `TarCA.beta` obtained from ggtree, which
#' `TarCA.pro` no longer stores. The result is identical to
#' `ape::node.depth.edgelength` with names `c(tip.label, node.label)`.
#'
#' @param Tree A tree of class "phylo" with node labels.
#'
#' @return A named numeric vector of root-to-node distances.
#' @noRd
fun.GetNodeDepth <- function(Tree){
  tmp.depth <- ape::node.depth.edgelength(Tree)
  names(tmp.depth) <- c(Tree$tip.label,Tree$node.label)
  tmp.depth
}
