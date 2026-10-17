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

Alternatively, install the dependency explicitly first:

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

### 2. Load demo data

Demo data can be downloaded from [LEAP_sourcedata](https://github.com/shadowdeng1994/LEAP_sourcedata). It includes:

- **Query data**: a reconstructed cell phylogeny (`Var.Tree`), tip annotations (`Var.Ann`), and a gene expression matrix (`Var.Expression`)
- **Reference data**: a reference expression matrix (`Data.Ref_Expression`) and metadata (`Data.Ref_metadata`)

```r
load("DemoData/Query.RData")
load("DemoData/Reference.RData")
```

### 3. Create TarCA object

Convert the phylogenetic tree and tip annotations into a TarCA object, the central data structure for all downstream analyses.

```r
Data.TarCA <- CreateTarCAObject(Var.Tree, Var.Ann)
```

### 4. Reconstruct internal node transcriptome

This is the core step. It automatically performs three sub-steps:

1. **LUG detection** — identifies lineage-specific upregulated genes using a permutation test (1000 shuffles by default)
2. **Phylogenetic depth conversion** — maps LUG activation events onto the phylogenetic depth of internal nodes
3. **Internal state inference** — reconstructs the binarized transcriptome of each internal node

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  Var.Expression,
  TopGene = 0.1,           # top fraction of cells binarized as "expressed"
  SeedFrom = 1,            # set (any integer) for reproducible LUG detection
  DirOut = "MyOutput/",    # output directory (created if missing)
  Threads = 90,
  SaveObject = TRUE
)
```

### 5. Binarize reference data

Convert the reference scRNA-seq dataset into a binarized matrix using the same
LUG genes, grouped by developmental time bins or other metadata categories.

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

### 6. Define transfer labels

Create combined labels from reference metadata (e.g., time bin and lineage) for
transferring to internal nodes.

```r
tmp.Label <- setNames(
  paste0(Data.Ref_metadata$embryo.time.bin, "|", Data.Ref_metadata$lineage),
  rownames(Data.Ref_metadata)
)
tmp.Label <- tmp.Label[colnames(Var.ReferenceMatrix)]   # align to the (grouped) reference columns
```

### 7. Run label transfer

Use Seurat's label transfer framework to assign developmental states to the reconstructed internal nodes.

```r
Var.Predictions <- RunLabelTransfer(
  Var.InternalMatrix,
  Var.ReferenceMatrix,
  tmp.Label,
  DirOut = "MyOutput/",
  SaveObject = TRUE,
  Prefix = "E7.5"          # optional; default: a random 4-character index
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
  SaveObject = TRUE,
  Prefix = "E7.5"          # use the same Prefix as RunLabelTransfer
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
| `fun.SetVerbose` | Set the log verbosity level (0 silent, 1 steps + timing, 2 details) |

## Output

All results are saved to the directory given by `DirOut` (default `"LEAPOutput/"`, created automatically) when `SaveObject = TRUE`:

- `Var.LUGs.rds` — detected lineage-specific upregulated genes with empirical P-values
- `Var.PhyloDepthOfLUGs.rds` — phylogenetic depth of LUG activation events
- `Var.InternalState.rds` — reconstructed binarized transcriptome of internal nodes (sparse `dgCMatrix`)
- `Var.RefMatrix.rds` — binarized reference matrix (sparse `dgCMatrix`)
- `Var.QuerySeurat.rds` / `Var.RefSeurat.rds` — Seurat objects for query and reference (label-independent; loaded automatically if already present)
- `Var.Anchors.rds` — transfer anchors (label-independent; loaded automatically if already present)
- `Var.Predictions_<index>.rds` — label transfer predictions, where `<index>` is `Prefix` or a random 4-character `[0-9A-Za-z]` string (always rewritten)
- `Var.PredictedState_<index>.rds` — final predicted states mapped onto the phylogeny, named the same way (always rewritten)
- `Var.Potency_<index>.rds` — per-internal-node potency (Shannon entropy of the clade's tip annotations), named the same way (always rewritten)

## Notes

### Using Seurat objects as input

Both the query and the reference can be supplied directly as **Seurat objects**
(their cells must match the tree tips). Extra cells are dropped; missing tips
raise an error. Choose the data with `Assay` and `Layer` (Seurat >= 5 uses a
`layer`; older Seurat uses the matching `slot`).

Query:

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  QuerySeurat,             # a Seurat object (cells aligned to Data.TarCA$Tree$tip.label)
  Assay = "RNA",           # default: DefaultAssay(QuerySeurat)
  Layer = "data",          # e.g. "data"/"counts"
  TopGene = 0.1,
  SeedFrom = 1,
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE
)
```

Reference (its metadata `obj[[]]` is used automatically, so `Metadata` can be omitted):

```r
RefSeurat <- subset(RefSeurat, features = intersect(rownames(Var.InternalMatrix), rownames(RefSeurat)))

Var.ReferenceMatrix <- GetReferenceMatrix(
  RefSeurat,                    # a Seurat object
  GroupBy = "embryo.time.bin",  # a column of RefSeurat's metadata
  Assay = "RNA",
  Layer = "data",
  TopGene = 0.1,
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE
)
```

### Running only the LUG step, or reusing a LUG set

Use `OnlyLUG = TRUE` to run **only** the LUG-detection step (it stops before
PhyloDepth conversion and internal-state inference) and return the detected LUGs:

```r
Var.LUGs <- GetInternalMatrix(
  Data.TarCA,
  Var.Expression,
  SeedFrom = 1,
  DirOut = "MyOutput/",
  Threads = 90,
  SaveObject = TRUE,
  OnlyLUG = TRUE          # return the LUG data.frame and stop
)
```

If you already have a set of LUG genes, pass them directly (a **character vector**
of gene symbols) to skip the expensive LUG-estimation step. The per-gene
`LUG/*.rds` files must already exist in `DirOut` (e.g. produced by a previous run):

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  DirOut = "MyOutput/",
  LUGs = Var.LUGs$Symbol[Var.LUGs$Sig]   # or any character vector of gene names
)
```

A supplied `LUGs` is an **input** (usually a subset of an already estimated
`Var.LUGs`), not a result, so it is **not** written to disk and never overwrites
`Var.LUGs.rds`. `Var.PhyloDepthOfLUGs.rds` and `Var.InternalState.rds` are still
written as usual.

### Re-running with a different label set

`Var.QuerySeurat.rds`, `Var.RefSeurat.rds` and `Var.Anchors.rds` depend only on
the two binarized matrices, **not** on the transfer labels, so `RunLabelTransfer`
loads them from `DirOut` when they already exist (set `Force = TRUE` to recompute
and overwrite them). Only the predictions depend on the labels.

The label-dependent outputs are **always rewritten**, and their names carry an
index: `<index>` is `Prefix` or — if `Prefix` is `NULL` (default) — a random
4-character `[0-9A-Za-z]` string.

| Function | Output |
|---|---|
| `RunLabelTransfer` | `Var.Predictions_<index>.rds` |
| `GetPredictedState` | `Var.PredictedState_<index>.rds` |
| `fun.CalculatePotency` | `Var.Potency_<index>.rds` |

Give the same `Prefix` to these calls so one label set produces one matching set
of files; a new `Prefix` writes new files, the same `Prefix` overwrites them.
Changing only the label vector therefore lets you run as many label sets as you
like without clobbering the others:

```r
Var.Predictions.E75 <- RunLabelTransfer(
  Var.InternalMatrix, Var.ReferenceMatrix, tmp.Label,
  DirOut = "MyOutput/", Prefix = "E7.5"
)
Var.PredictedState.E75 <- GetPredictedState(
  Data.TarCA, Var.Predictions.E75,
  DirOut = "MyOutput/", Prefix = "E7.5"
)
Var.Potency.E75 <- fun.CalculatePotency(
  Data.TarCA, Var.PredictedState.E75$NodeLabel,
  DirOut = "MyOutput/", Prefix = "E7.5"
)
```

With the default random index (no `Prefix`) every run writes a new file, so pass
`Prefix` if you want a stable, recognisable name.

### Parallelism

Pass `Threads` to the functions above; each call runs its gene- and
permutation-level work on a `future` **multisession** backend with that many
workers (works on Windows, macOS and Linux, and avoids the fork-related failures
of `parallel::mclapply` inside RStudio). Alternatively, set your own plan once,
which then takes precedence over `Threads`:

```r
future::plan(future::multisession, workers = 8)
```

LUG estimation sends the expression data to the workers in chunks (the data
rides inside `future`'s `X`, so each worker only receives its own chunk instead
of a copy of the whole matrix). The default is **one chunk per worker**
(`ceiling(nGenes / Threads)`), which is optimal when the per-gene cost is
uniform. Control it with `ChunkSize = <genes per chunk>` if needed — smaller
chunks (with dynamic scheduling) only help when per-gene cost varies a lot.

### Large datasets (`future.globals.maxSize`)

With `multisession`, the objects a worker needs are *copied* to the workers, so a
large query/reference matrix can exceed `future`'s default 500 MiB limit and raise
`Error ... The total size of the ... globals ... exceeds the maximum allowed size 500.00 MiB`.
Estimate the largest object that will be exported and raise the limit before
calling the analysis functions:

```r
# size of the biggest object sent to the workers (e.g. the expression matrix)
format(object.size(Var.Expression), units = "MiB")
# a dense genes x cells matrix is about nrow(Var.Expression) * ncol(Var.Expression) * 8 bytes;
# for a sparse (dgCMatrix) input use object.size() directly

options(future.globals.maxSize = 2 * 1024^3)   # set to ~1.5-2x the largest exported object
# only if you are certain the transfer is intended:
options(future.globals.maxSize = Inf)
```

For very large data on Linux/macOS, `future::plan(future::multicore)` shares
memory between workers instead of copying it, but `multicore` forks and can hit
the same RStudio issues as `parallel::mclapply`.

### Verbosity and timing

Progress output is controlled by a verbosity level (`options(LEAP.verbose)`,
default `1`). Every step reports its elapsed time so you can see where time is
spent. Set it **before** running the analysis:

```r
fun.SetVerbose(2)              # 0 = silent, 1 = main steps + timing (default), 2 = per-step details
# or, equivalently:
options(LEAP.verbose = 2)
```

For example, `GetInternalMatrix` prints the time of each sub-step, and at level
`2` the per-`Total` shuffling times:

```
[..] ==> Extracting LUG.
[..] --> Running LUG estimator.
[..] --> LUG estimator done in  0.35 s.
[..] --> Shuffling ( 10 times).
[..] --> Shuffling done in  0.09 s.
[..] --> LUG detection done in  0.22 s. ( 107  significant)
[..] ==> Converting to PhyloDepth done in  0.32 s.
[..] ==> Inferring internal state done in  0.03 s.
[..] ==> GetInternalMatrix total  1.22 s.
```

### Reproducibility, thresholds and caching

- **`SeedFrom`** makes the permutation test (empirical P-values) reproducible;
  leave it `NULL` (default) to use random seeds.
- **`TopGene`** is the top fraction of cells (within each group) treated as
  "highly expressed" when binarizing. Use the **same value for the query and the
  reference** so the two binarized matrices are comparable (default `0.1`;
  smaller is stricter, larger is more permissive).
- **`DirOut`** is where results and the cached `LUG/` and `Shuffle/` files are
  written. Use the same `DirOut` across steps so the caches can be reused.
- Already-computed `LUG/*.rds` and `Shuffle/*.rds` are **skipped by default**;
  set **`Force = TRUE`** to recompute them (e.g. after changing `TopGene` or
  `SeedFrom`). `Redo` is kept as an alias for `Force`.
- **`RunLabelTransfer`** likewise reuses `Var.QuerySeurat.rds`,
  `Var.RefSeurat.rds` and `Var.Anchors.rds` from `DirOut` when they exist
  (`Force = TRUE` recomputes them). Its label-dependent output,
  `Var.Predictions_<Prefix|random>.rds`, is always rewritten, as is
  `GetPredictedState`'s `Var.PredictedState_<Prefix|random>.rds`.

## Contributing
### Contributors
Shanjun Deng, shadowdeng1994@gmail.com.

## Citations
When using TarCA please cite:

We're on track to publish something really exciting! **Star the repo and stay tuned!​** :smile:
