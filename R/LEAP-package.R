utils::globalVariables(c(
  "predicted.id","TransLabel","AncLabel","isTip","isRoot","Total",".",
  "Np","EmpiricalP","TipLabel","Nor_x","Symbol","Node","Node_x","Symbol_x","isAct",
  "Sig","TipAnn","N",".N"
))

#' LEAP.pro: LUG-Encoded Ancestral Projection
#'
#' A phylogeny-based computational framework for inferring the transcriptional
#' states of internal nodes from cell phylogenies via lineage-specific
#' upregulated genes (LUGs).
#'
#' @importFrom TarCA.pro Np_Estimator
#' @importFrom ape node.depth.edgelength
#' @importFrom matrixStats rowRanks
#' @importFrom data.table as.data.table
#' @importFrom dplyr filter mutate select arrange bind_rows left_join group_by percent_rank bind_cols
#' @importFrom tidyr spread gather replace_na
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom Seurat CreateSeuratObject FindTransferAnchors TransferData NormalizeData ScaleData RunPCA RunUMAP
#' @importFrom stats median setNames p.adjust pchisq
#' @importFrom utils flush.console
"_PACKAGE"
