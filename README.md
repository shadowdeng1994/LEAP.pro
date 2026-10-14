# LUG-Encoded Ancestral Projection (LEAP)

**Package:** `LEAP.pro` — an R implementation of LEAP (maintained separately from the `LEAP` package).

A phylogeny-based computational framework for inferring the transcriptional states of unobserved internal nodes from cell phylogenies via lineage-specific upregulated genes (LUGs).

## Overview

LEAP integrates cell lineage tracing with single-cell RNA sequencing (scRNA-seq) to reconstruct the transcriptomes of ancestral internal nodes in a developmental cell lineage tree. It identifies a set of lineage-committed upregulated genes (LUGs) from terminal cell transcriptomes, maps them back onto internal nodes of the phylogeny, and enables lineage-informed longitudinal analysis of cell state dynamics throughout development.

The method facilitates the discovery of subtle fate specializations and hidden developmental states that lack prominent transcriptional differentiation, supporting the construction of lineage-resolved cell atlases in complex organisms.

<img width="1924" height="849" alt="image" src="https://github.com/user-attachments/assets/f922c00b-b1d6-4bf4-b22e-2dd76ca3b3a2" />

## Installation

The non-CRAN dependency [TarCA.pro](https://github.com/shadowdeng1994/TarCA.pro)
is resolved automatically via the `Remotes` field, so `LEAP.pro` can be installed
in a single step with [pak](https://pak.r-lib.org/):

```r
# install.packages("pak")
pak::pkg_install("shadowdeng1994/LEAP.pro")
```

Alternatively, install the dependency explicitly before `LEAP.pro`:

```r
pak::pkg_install("shadowdeng1994/TarCA.pro")
pak::pkg_install("shadowdeng1994/LEAP.pro")
```

## Dependencies

- R (>= 4.0.0)
- [TarCA.pro](https://github.com/shadowdeng1994/TarCA.pro) (tree/ExTree construction)
- dplyr, tidyr, tibble, magrittr, data.table
- future, future.apply (parallel backend), matrixStats, Matrix, ape
- Seurat

## Workflow

### 1. Load packages

```r
library(LEAP.pro)
library(TarCA.pro)
library(future)
library(Seurat)
```

**Parallelism.** Pass `Threads` to the functions below; each call runs its gene-
and permutation-level work on a `future` **multisession** backend with that many
workers (works on Windows, macOS and Linux, and avoids the fork-related failures
of `parallel::mclapply` inside RStudio). Alternatively, set your own plan once,
which then takes precedence over `Threads`:

```r
future::plan(future::multisession, workers = 8)
```

**Large datasets (`future.globals.maxSize`).** With `multisession`, the objects a
worker needs are *copied* to the workers, so a large query/reference matrix can
exceed `future`'s default 500 MiB limit and raise
`Error ... The total size of the ... globals ... exceeds the maximum allowed size 500.00 MiB`.
Estimate the largest object that will be exported and raise the limit before
calling the analysis functions:

```r
# size of the biggest object sent to the workers (e.g. the expression matrix)
format(object.size(Var.Expression), units = "MiB")
# a dense genes x cells matrix is about nrow(Var.Expression) * ncol(Var.Expression) * 8 bytes;
# for a sparse (dgCMatrix) input use object.size() directly, or vapply over a few genes

options(future.globals.maxSize = 2 * 1024^3)   # 2 GiB: set to ~1.5-2x the largest exported object
# only if you are certain the transfer is intended:
options(future.globals.maxSize = Inf)
```

The option must be set **before** calling `GetInternalMatrix()` /
`GetReferenceMatrix()`. For very large data on Linux/macOS, `future::plan(future::multicore)`
shares memory between workers instead of copying it, but `multicore` forks and can
hit the same RStudio issues as `parallel::mclapply`.

### 2. Load demo data

Demo data can be downloaded from [LEAP_sourcedata](https://github.com/shadowdeng1994/LEAP_sourcedata). It includes:

- **Query data**: a reconstructed cell phylogeny (`Var.Tree`), tip annotations (`Var.Ann`), and a gene expression matrix (`Var.Expression`)
- **Reference data**: a reference expression matrix (`Data.Ref_Expression`) and metadata (`Data.Ref_metadata`) for label transfer

```r
load("DemoData/Query.RData")
load("DemoData/Reference.RData")
```

### 3. Create TarCA object

Convert the phylogenetic tree and tip annotations into a TarCA object, which serves as the central data structure for all downstream analyses.

```r
Data.TarCA <- CreateTarCAObject(Var.Tree, Var.Ann)
```

### 4. Reconstruct internal node transcriptome

This is the core step of LEAP. It performs three sub-steps automatically:

1. **LUG detection** — identifies lineage-specific upregulated genes using a permutation test (1000 shuffles by default)
2. **Phylogenetic depth conversion** — maps LUG activation events onto the phylogenetic depth of internal nodes
3. **Internal state inference** — reconstructs the binarized transcriptome of each internal node

The permutation test uses `SeedFrom` to make the empirical P-values reproducible; leave it `NULL` (default) to use random seeds. `TopGene` is the top fraction of cells (within each group) treated as "highly expressed" when binarizing; use the **same value for the query and the reference** so the two binarized matrices are comparable (default `0.1`; smaller is stricter, larger is more permissive). Results are written to `DirOut` (default `"LEAPOutput/"`, created automatically); use the same `DirOut` across steps so the cached `LUG/` and `Shuffle/` files can be reused.

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  Var.Expression,
  TopGene = 0.1,           # top fraction of cells binarized as "expressed" (default 0.1)
  SeedFrom = 1,            # set (any integer) for reproducible LUG detection
  DirOut = "MyOutput/",    # output directory (created if missing)
  Threads = 90,
  SaveObject = TRUE
)
```

`ExpreMatrix` can also be a **Seurat object** whose cells match the tree tips (extra cells are dropped; missing tips raise an error). Pick the assay and layer/slot with `Assay` and `Layer`:

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  QuerySeurat,             # a Seurat object (cells aligned to Data.TarCA$Tree$tip.label)
  Assay = "RNA",           # assay to use (default: DefaultAssay(QuerySeurat))
  Layer = "data",          # Seurat >= 5 layer, e.g. "data"/"counts"; older Seurat uses the matching slot
  TopGene = 0.1,
  SeedFrom = 1,
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE
)
```

### 5. Binarize reference data

Convert the reference scRNA-seq dataset into a binarized matrix using the same LUG genes, grouped by developmental time bins or other metadata categories.

```r
tmp.overlapLUG <- intersect(rownames(Var.InternalMatrix), rownames(Data.Ref_Expression))

Var.ReferenceMatrix <- GetReferenceMatrix(
  Data.Ref_Expression[tmp.overlapLUG, ],
  Data.Ref_metadata,
  GroupBy = "embryo.time.bin",
  TopGene = 0.1,           # must match the query's TopGene
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE
)
```

`ExpreMatrix` can instead be a **Seurat object**; its metadata (`obj[[]]`) is then used automatically and `Metadata` can be omitted. Subset to the same LUG genes first, and choose the data with `Assay`/`Layer`:

```r
RefSeurat <- subset(RefSeurat, features = intersect(rownames(Var.InternalMatrix), rownames(RefSeurat)))

Var.ReferenceMatrix <- GetReferenceMatrix(
  RefSeurat,               # a Seurat object
  GroupBy = "embryo.time.bin",  # a column of RefSeurat's metadata
  Assay = "RNA",           # assay to use (default: DefaultAssay(RefSeurat))
  Layer = "data",          # Seurat >= 5 layer, e.g. "data"/"counts"; older Seurat uses the matching slot
  TopGene = 0.1,
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE
)
```

### 6. Define transfer labels

Create combined labels from reference metadata (e.g., time bin and lineage) for transferring to internal nodes.

```r
tmp.Label <- setNames(
  paste0(Data.Ref_metadata$embryo.time.bin, "|", Data.Ref_metadata$lineage),
  rownames(Data.Ref_metadata)
)
tmp.Label <- tmp.Label[colnames(Var.ReferenceMatrix)]
```

### 7. Run label transfer

Use Seurat's label transfer framework to assign developmental states to the reconstructed internal nodes.

```r
Var.Predictions <- RunLabelTransfer(
  Var.InternalMatrix,
  Var.ReferenceMatrix,
  tmp.Label,
  DirOut = "MyOutput/",
  SaveObject = TRUE
)
```

### 8. Extract predicted states

Map the transferred labels onto the cell phylogeny, propagating annotations to unsolved descendant nodes.

```r
Var.PredictedState <- GetPredictedState(
  Data.TarCA,
  Var.Predictions,
  RemoveRoot = "Root",
  DirOut = "MyOutput/",
  SaveObject = TRUE
)
```

## Key Functions

| Function | Description |
|---|---|
| `CreateTarCAObject` | Create a TarCA object from a phylogenetic tree and tip annotations |
| `GetInternalMatrix` | Reconstruct LUG-transcriptome of internal nodes (main workflow) |
| `GetReferenceMatrix` | Binarize reference dataset by expression percentiles within groups |
| `RunLabelTransfer` | Transfer labels from reference to query data via Seurat |
| `GetPredictedState` | Map inferred states to nodes and their descendants |
| `fun.ExtractLUG` | Infer LUGs based on permutation test |
| `fun.CalculatePotency` | Calculate Shannon entropy-based potency for each internal node |
| `fun.AssignLaterState` | Propagate predicted states through the internal node hierarchy |

## Output

All results are saved to the directory given by `DirOut` (default `"LEAPOutput/"`, created automatically) when `SaveObject = TRUE`:

- `Var.LUGs.rds` — detected lineage-specific upregulated genes with empirical P-values
- `Var.PhyloDepthOfLUGs.rds` — phylogenetic depth of LUG activation events
- `Var.InternalState.rds` — reconstructed binarized transcriptome of internal nodes
- `Var.RefMatrix.rds` — binarized reference matrix
- `Var.QuerySeurat.rds` / `Var.RefSeurat.rds` — Seurat objects for query and reference
- `Var.Anchors.rds` — transfer anchors
- `Var.Predictions.rds` — label transfer predictions
- `Var.PredictedState.rds` — final predicted states mapped onto the phylogeny

## Contributing
### Contributors
Shanjun Deng, shadowdeng1994@gmail.com.

## Citations
When using TarCA please cite:

We're on track to publish something really exciting! **Star the repo and stay tuned!​** :smile:
