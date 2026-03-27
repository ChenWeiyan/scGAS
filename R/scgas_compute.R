#' Compute Single-Cell Gene Activation Scores via Network Propagation
#'
#' Takes the Metacell-level GAS values produced by
#' \code{\link{scgas_train_models}} and propagates them to single-cell
#' resolution in two steps:
#' \enumerate{
#'   \item **Initialisation (scGAS0).**  Each non-centroid cell inherits a
#'     weighted average of the GAS values from its three nearest Metacell
#'     centroids.  Weights are computed with a Gaussian kernel on Euclidean
#'     distances in the LSI embedding space.
#'   \item **Network propagation (scGAS).**  A mutual-nearest-neighbour (MNN)
#'     graph is built over the LSI embedding.  The top \code{seed_fraction}
#'     cells (those with the highest scGAS0) are used as seed nodes for a
#'     random-walk propagation via \code{SCAVENGE::randomWalk_sparse}.  The
#'     result is outlier-clipped and min-max scaled to \eqn{[0, 1]}.
#' }
#'
#' The function returns the scGAS matrix and, optionally, adds it as a new
#' \code{scGAS} assay inside the supplied \code{SeuratObject}.
#'
#' @param seurat_obj \code{SeuratObject} returned by
#'   \code{\link{scgas_preprocess}}.  Must contain an \code{lsi} reduction and
#'   \code{seurat_clusters} metadata.
#' @param trained_models Named list returned by
#'   \code{\link{scgas_train_models}}.  The attributes
#'   \code{"ATAC_CL_lognorm"} and \code{"cl_levels"} must be present.
#' @param lsi_dims Integer vector of LSI dimensions used for the embedding.
#'   Must match those used in \code{\link{scgas_preprocess}}.  Default
#'   \code{2:30}.
#' @param knn_k Integer.  Number of mutual nearest neighbours.  Default
#'   \code{30}.
#' @param n_metacell_refs Integer.  Number of nearest Metacell centroids used
#'   to initialise each cell's scGAS0.  Default \code{3}.
#' @param seed_fraction Numeric in (0, 1].  Fraction of cells (ranked by
#'   scGAS0) to use as seed nodes in the propagation.  Default \code{0.10}.
#' @param gamma Numeric.  Restart probability for the random walk.  Default
#'   \code{0.05}.
#' @param outlier_quantile Numeric.  Upper quantile used to cap outlier
#'   propagation scores before min-max scaling.  Default \code{0.95}.
#' @param n_cores Integer.  Cores for \code{parallel::mclapply}.  Default
#'   \code{1}.
#' @param verbose Logical.  Default \code{TRUE}.
#' @param metacell_membership Named integer vector of Metacell assignments
#'   as returned by \code{\link{scgas_metacell}}.  When \code{NULL} (default),
#'   falls back to \code{seurat_obj$seurat_clusters} for backward compatibility.
#'
#' @return A genes x cells numeric matrix of scGAS values scaled to
#'   \eqn{[0, 1]}.  Row names are gene symbols; column names are cell
#'   barcodes.
#'
#' @seealso \code{\link{scgas_add_assay}} to embed the matrix into a
#'   \code{SeuratObject}.
#'
#' @examples
#' \dontrun{
#' scgas_mat <- scgas_compute(
#'   seurat_obj     = obj,
#'   trained_models = models,
#'   lsi_dims       = 2:30,
#'   knn_k          = 30,
#'   seed_fraction  = 0.10,
#'   n_cores        = 8
#' )
#' }
#'
#' @importFrom MatrixGenerics colMeans2
#' @importFrom parallel mclapply
#' @importFrom SCAVENGE getmutualknn randomWalk_sparse capOutlierQuantile
#'   max_min_scale
#' @importFrom stats dist
#' @export
scgas_compute <- function(seurat_obj,
                          trained_models,
                          lsi_dims         = 2:30,
                          knn_k            = 30L,
                          n_metacell_refs  = 3L,
                          seed_fraction    = 0.10,
                          gamma            = 0.05,
                          outlier_quantile = 0.95,
                          n_cores          = 1L,
                          verbose          = TRUE,
                          metacell_membership) {

  ## ── Validation ─────────────────────────────────────────────────────────────
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    is.list(trained_models), !is.null(names(trained_models)),
    !is.null(attr(trained_models, "ATAC_CL_lognorm")),
    is.numeric(seed_fraction), seed_fraction > 0, seed_fraction <= 1,
    is.numeric(gamma), gamma > 0
  )
  if (!"lsi" %in% names(seurat_obj@reductions))
    stop("[scGAS] 'lsi' reduction not found. Run scgas_preprocess() first.")

  ## ── Extract LSI embedding ──────────────────────────────────────────────────
  embeddings <- seurat_obj@reductions$lsi@cell.embeddings[, lsi_dims,
                                                           drop = FALSE]
  n_cells    <- nrow(embeddings)
  if (!is.null(metacell_membership)) {
    cl_labels <- metacell_membership[colnames(seurat_obj)]
  } else {
    cl_labels <- seurat_obj$seurat_clusters
  }
  cl_levels  <- attr(trained_models, "cl_levels")
  sel_genes  <- names(trained_models)

  if (verbose) message("[scGAS] Building mutual-KNN graph (k = ", knn_k, ")")
  mutualknn <- SCAVENGE::getmutualknn(embeddings, knn_k)

  ## ── Identify Metacell centroid cells ──────────────────────────────────────
  if (verbose) message("[scGAS] Identifying Metacell centroids")
  centroid_cells <- .find_centroids(embeddings, cl_labels, cl_levels)
  ## centroid_cells: named character vector [Cluster-X] -> barcode

  ## ── Pre-compute Gaussian-kernel distances: all cells to centroids ──────────
  if (verbose) message("[scGAS] Computing cell-to-centroid distances")
  other_cells   <- base::setdiff(rownames(embeddings), centroid_cells)
  dist_all      <- as.matrix(stats::dist(embeddings,
                                          method = "euclidean",
                                          diag   = TRUE,
                                          upper  = TRUE))
  ## Bandwidth: distance to knn_k-th nearest neighbour
  alpha_vec     <- apply(dist_all, 1,
                         function(x) sort(x, decreasing = FALSE)[knn_k])

  ## Build the cell -> top-3-centroids reference
  cells_ref <- .build_cell_refs(
    other_cells   = other_cells,
    centroid_cells = centroid_cells,
    dist_all      = dist_all,
    alpha_vec     = alpha_vec,
    n_metacell_refs = n_metacell_refs
  )

  ## ── Per-gene scGAS computation ────────────────────────────────────────────
  if (verbose) message("[scGAS] Computing scGAS for ", length(sel_genes),
                       " genes")

  ## Chunk genes for memory-friendly parallel execution
  chunk_size  <- 500L
  n_genes     <- length(sel_genes)
  n_chunks    <- max(1L, round(n_genes / chunk_size))
  gene_chunks <- split(seq_len(n_genes),
                       sample(factor(seq_len(n_genes) %% n_chunks)))

  results     <- vector("list", n_chunks)
  for (chunk_i in seq_len(n_chunks)) {
    idx_v      <- gene_chunks[[chunk_i]]
    chunk_genes <- sel_genes[idx_v]

    chunk_res <- parallel::mclapply(
      X   = as.list(chunk_genes),
      FUN = .compute_one_gene,
      trained_models  = trained_models,
      centroid_cells  = centroid_cells,
      other_cells     = other_cells,
      cells_ref       = cells_ref,
      n_cells         = n_cells,
      cell_names      = colnames(seurat_obj),
      mutualknn       = mutualknn,
      seed_fraction   = seed_fraction,
      gamma           = gamma,
      outlier_quantile = outlier_quantile,
      mc.cores        = n_cores
    )
    names(chunk_res)   <- chunk_genes
    results[[chunk_i]] <- chunk_res
    gc()
  }

  ## ── Assemble output matrix ─────────────────────────────────────────────────
  if (verbose) message("[scGAS] Assembling result matrix")
  scgas_mat  <- matrix(
    NA_real_,
    nrow     = n_genes,
    ncol     = ncol(seurat_obj),
    dimnames = list(sel_genes, colnames(seurat_obj))
  )
  for (chunk_i in seq_len(n_chunks)) {
    for (g in names(results[[chunk_i]])) {
      res_g <- results[[chunk_i]][[g]]
      if (!is.null(res_g))
        scgas_mat[g, ] <- res_g$scGAS
    }
  }

  if (verbose) message("[scGAS] Done. scGAS matrix: ",
                       nrow(scgas_mat), " genes x ",
                       ncol(scgas_mat), " cells.")
  return(scgas_mat)
}


#' Add scGAS Matrix as an Assay Inside a SeuratObject
#'
#' Convenience wrapper that stores the genes x cells scGAS matrix as a new
#' \code{"scGAS"} assay and, optionally, re-runs PCA and UMAP on the scGAS
#' values to produce a cell embedding.
#'
#' @param seurat_obj \code{SeuratObject}.
#' @param scgas_mat Genes x cells numeric matrix as returned by
#'   \code{\link{scgas_compute}}.
#' @param run_dim_reduction Logical.  If \code{TRUE} (default) run
#'   \code{FindVariableFeatures}, \code{ScaleData}, \code{RunPCA}, and
#'   \code{RunUMAP} on the new assay.
#' @param n_pcs Integer.  Number of PCs for PCA and UMAP.  Default \code{30}.
#' @param n_neighbors Integer.  \code{n.neighbors} passed to
#'   \code{RunUMAP}.  Default \code{30}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return The \code{SeuratObject} with a new \code{scGAS} assay (and
#'   optionally a \code{scgasumap} reduction).
#'
#' @examples
#' \dontrun{
#' obj <- scgas_add_assay(obj, scgas_mat)
#' }
#'
#' @importFrom Seurat CreateAssayObject DefaultAssay FindVariableFeatures
#'   ScaleData RunPCA RunUMAP
#' @export
scgas_add_assay <- function(seurat_obj,
                            scgas_mat,
                            run_dim_reduction = TRUE,
                            n_pcs             = 30L,
                            n_neighbors       = 30L,
                            verbose           = TRUE) {
  stopifnot(inherits(seurat_obj, "Seurat"), is.matrix(scgas_mat))

  seurat_obj[["scGAS"]] <- Seurat::CreateAssayObject(counts = scgas_mat)
  seurat_obj@assays$scGAS@data <- scgas_mat

  if (run_dim_reduction) {
    if (verbose) message("[scGAS] Running dimensionality reduction on scGAS")
    Seurat::DefaultAssay(seurat_obj) <- "scGAS"
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = verbose)
    seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = verbose)
    seurat_obj <- Seurat::RunPCA(
      seurat_obj,
      features    = Seurat::VariableFeatures(seurat_obj),
      reduction.name = "scgaspca",
      verbose     = verbose
    )
    seurat_obj <- Seurat::RunUMAP(
      seurat_obj,
      reduction      = "scgaspca",
      dims           = seq_len(n_pcs),
      n.components   = 2L,
      n.neighbors    = n_neighbors,
      reduction.name = "scgasumap",
      verbose        = verbose
    )
  }
  return(seurat_obj)
}


## ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.find_centroids <- function(embeddings, cl_labels, cl_levels) {
  centroids <- character(length(cl_levels))
  names(centroids) <- paste0("Cluster-", cl_levels)
  for (cl in cl_levels) {
    sel_cells   <- which(cl_labels == cl)
    cl_emb      <- embeddings[sel_cells, , drop = FALSE]
    cl_mean     <- MatrixGenerics::colMeans2(cl_emb)
    dists       <- apply(cl_emb, 1,
                         function(x) sqrt(sum((x - cl_mean)^2)))
    centroids[paste0("Cluster-", cl)] <- names(sort(dists))[1]
  }
  return(centroids)
}

#' @keywords internal
.build_cell_refs <- function(other_cells,
                             centroid_cells,
                             dist_all,
                             alpha_vec,
                             n_metacell_refs) {
  refs <- vector("list", length(other_cells))
  names(refs) <- other_cells

  for (cell in other_cells) {
    alpha     <- alpha_vec[cell]
    raw_dists <- dist_all[cell, centroid_cells]
    weights   <- exp(-(raw_dists / alpha)^2)
    weights   <- weights / sum(weights)
    top_idx   <- order(weights, decreasing = TRUE)[
      seq_len(min(n_metacell_refs, length(weights)))]
    top_wgt   <- weights[top_idx]
    refs[[cell]] <- top_wgt / sum(top_wgt)  # re-normalise
  }
  return(refs)
}

#' @keywords internal
.compute_one_gene <- function(gene_name,
                              trained_models,
                              centroid_cells,
                              other_cells,
                              cells_ref,
                              n_cells,
                              cell_names,
                              mutualknn,
                              seed_fraction,
                              gamma,
                              outlier_quantile) {
  tryCatch({
    model_g <- trained_models[[gene_name]]
    gas_mc  <- model_g$GAS
    names(gas_mc) <- unname(centroid_cells)  # centroid barcode -> GAS

    ## Initialise scGAS0
    scgas0 <- stats::setNames(numeric(length(cell_names)), cell_names)
    scgas0[unname(centroid_cells)] <- gas_mc

    ## Weighted average from nearest centroids for non-centroid cells
    for (cell in other_cells) {
      ref_wgt       <- cells_ref[[cell]]
      scgas0[cell]  <- sum(ref_wgt * gas_mc[names(ref_wgt)])
    }

    ## Network propagation
    seed_n    <- max(1L, floor(seed_fraction * length(scgas0)))
    seed_idx  <- rownames(mutualknn)[
      order(scgas0, decreasing = TRUE)[seq_len(seed_n)]
    ]
    np_score  <- SCAVENGE::randomWalk_sparse(
      intM      = mutualknn,
      queryCells = seed_idx,
      gamma     = gamma
    )
    np_score2 <- SCAVENGE::capOutlierQuantile(np_score, outlier_quantile)
    scgas     <- as.numeric(SCAVENGE::max_min_scale(np_score2))
    names(scgas) <- rownames(mutualknn)
    ## Ensure original cell order
    scgas <- scgas[cell_names]

    list(scGAS0 = scgas0, scGAS = scgas)
  }, error = function(e) NULL)
}
