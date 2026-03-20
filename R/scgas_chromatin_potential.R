#' Build and Visualise the Chromatin Potential Field
#'
#' The Chromatin Potential Field (CPF) links chromatin-level activation (scGAS)
#' to transcriptomic state.  For every cell \eqn{X}, the \eqn{K} cells whose
#' RNA expression profile is most correlated with \eqn{X}'s scGAS vector are
#' identified.  The centroid of those \eqn{K} neighbours in the low-dimensional
#' embedding defines an arrow whose direction and magnitude represent the
#' predicted direction of differentiation.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Select high-variable genes that are informative at both the scGAS
#'     and RNA levels.
#'   \item For each cell \eqn{X}, compute Spearman correlation between
#'     \eqn{X}'s scGAS profile and every other cell's RNA profile over the
#'     selected gene set.
#'   \item Retain the top \code{knn_k} correlated cells as the CPF neighbours
#'     of \eqn{X}.
#'   \item Project arrows onto the 2-D embedding supplied in
#'     \code{embedding_coords}.
#' }
#'
#' @param seurat_obj \code{SeuratObject} containing an \code{RNA} assay and
#'   a \code{scGAS} assay (added by \code{\link{scgas_add_assay}}).
#' @param scgas_mat Genes x cells scGAS matrix.  If \code{NULL} (default) the
#'   function attempts to read the \code{scGAS} assay from \code{seurat_obj}.
#' @param rna_mat Genes x cells RNA log-normalised matrix.  If \code{NULL}
#'   the function uses \code{seurat_obj@assays$RNA@data}.
#' @param embedding_coords Numeric matrix (cells x 2) of 2-D coordinates used
#'   to draw arrows.  Row names must match cell barcodes.  If \code{NULL} the
#'   function extracts the first two components of the \code{umap} or
#'   \code{scgasumap} reduction.
#' @param n_hvg Integer.  Number of highly variable genes to use.  Default
#'   \code{2000}.
#' @param knn_k Integer.  Number of RNA neighbours used to define each cell's
#'   CPF.  Default \code{10}.
#' @param n_cores Integer.  Default \code{1}.
#' @param arrow_scale Numeric.  Multiplicative scaling applied to arrow lengths
#'   for visualisation.  Default \code{1}.
#' @param point_colour Character.  Column name in \code{seurat_obj@meta.data}
#'   used to colour cells in the plot.  Default \code{"seurat_clusters"}.
#' @param point_size Numeric.  Size of background cell points.  Default
#'   \code{0.5}.
#' @param arrow_alpha Numeric in (0, 1].  Transparency of arrows.
#'   Default \code{0.7}.
#' @param title Character.  Plot title.  Default \code{"Chromatin Potential Field"}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{cpf_df}}{Data frame with one row per cell containing
#'       columns \code{x}, \code{y} (start coordinates), \code{dx}, \code{dy}
#'       (arrow components), and \code{magnitude}.}
#'     \item{\code{plot}}{A \code{ggplot2} object showing cells as points
#'       and CPF as arrows.}
#'     \item{\code{hvg}}{Character vector of genes used.}
#'   }
#'
#' @examples
#' \dontrun{
#' cpf <- scgas_chromatin_potential(
#'   seurat_obj    = obj,
#'   n_hvg         = 2000,
#'   knn_k         = 10,
#'   point_colour  = "celltype"
#' )
#' print(cpf$plot)
#' }
#'
#' @importFrom Seurat VariableFeatures FindVariableFeatures
#' @importFrom Matrix t
#' @importFrom stats cor setNames
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_colour_discrete
#'   theme_classic labs theme element_text
#' @export
scgas_chromatin_potential <- function(seurat_obj,
                                     scgas_mat        = NULL,
                                     rna_mat          = NULL,
                                     embedding_coords = NULL,
                                     n_hvg            = 2000L,
                                     knn_k            = 10L,
                                     n_cores          = 1L,
                                     arrow_scale      = 1,
                                     point_colour     = "seurat_clusters",
                                     point_size       = 0.5,
                                     arrow_alpha      = 0.7,
                                     title            = "Chromatin Potential Field",
                                     verbose          = TRUE) {

  ## ── Gather matrices ────────────────────────────────────────────────────────
  if (is.null(scgas_mat)) {
    if (!"scGAS" %in% names(seurat_obj@assays))
      stop("[scGAS] scgas_mat is NULL and no 'scGAS' assay found in ",
           "seurat_obj. Run scgas_add_assay() first.")
    scgas_mat <- seurat_obj@assays$scGAS@data
    if (!is.matrix(scgas_mat))
      scgas_mat <- as.matrix(scgas_mat)
  }

  if (is.null(rna_mat)) {
    if (!"RNA" %in% names(seurat_obj@assays))
      stop("[scGAS] No RNA assay found in seurat_obj.")
    rna_mat <- as.matrix(seurat_obj@assays$RNA@data)
  }

  ## ── 2-D embedding ──────────────────────────────────────────────────────────
  if (is.null(embedding_coords)) {
    red_names <- names(seurat_obj@reductions)
    if ("scgasumap" %in% red_names) {
      embedding_coords <-
        seurat_obj@reductions$scgasumap@cell.embeddings[, 1:2]
    } else if ("umap" %in% red_names) {
      embedding_coords <-
        seurat_obj@reductions$umap@cell.embeddings[, 1:2]
    } else {
      stop("[scGAS] No UMAP reduction found. Supply embedding_coords.")
    }
  }
  colnames(embedding_coords) <- c("UMAP_1", "UMAP_2")

  ## ── Select highly variable genes common to both matrices ──────────────────
  if (verbose) message("[scGAS] Selecting HVGs for CPF computation")
  common_genes <- intersect(rownames(scgas_mat), rownames(rna_mat))
  if (length(common_genes) == 0)
    stop("[scGAS] No genes shared between scgas_mat and rna_mat.")

  ## Variable genes based on scGAS variance
  scgas_vars    <- apply(scgas_mat[common_genes, , drop = FALSE], 1, stats::var)
  rna_vars      <- apply(rna_mat[common_genes, , drop = FALSE], 1, stats::var)
  hvg_scgas     <- names(sort(scgas_vars, decreasing = TRUE))[
    seq_len(min(n_hvg, length(common_genes)))]
  hvg_rna       <- names(sort(rna_vars, decreasing = TRUE))[
    seq_len(min(n_hvg, length(common_genes)))]
  hvg           <- intersect(hvg_scgas, hvg_rna)

  if (length(hvg) < 10)
    warning("[scGAS] Only ", length(hvg), " HVGs shared between scGAS and ",
            "RNA. Results may be noisy.")
  if (verbose) message("[scGAS] Using ", length(hvg), " HVGs")

  ## Subset matrices
  gas_hvg  <- t(scgas_mat[hvg, , drop = FALSE])   # cells x hvg
  rna_hvg  <- t(rna_mat[hvg, , drop = FALSE])      # cells x hvg
  all_cells <- colnames(seurat_obj)

  ## ── Compute CPF arrows in parallel ────────────────────────────────────────
  if (verbose) message("[scGAS] Computing CPF arrows (k = ", knn_k, ")")
  cpf_list <- parallel::mclapply(
    X        = as.list(all_cells),
    FUN      = .cpf_one_cell,
    gas_hvg  = gas_hvg,
    rna_hvg  = rna_hvg,
    embedding = embedding_coords,
    knn_k    = knn_k,
    mc.cores = n_cores
  )
  names(cpf_list) <- all_cells

  ## Assemble result data frame
  cpf_df <- do.call(rbind, lapply(all_cells, function(cell) {
    r <- cpf_list[[cell]]
    if (is.null(r)) return(NULL)
    data.frame(
      cell      = cell,
      x         = embedding_coords[cell, 1],
      y         = embedding_coords[cell, 2],
      dx        = r$dx * arrow_scale,
      dy        = r$dy * arrow_scale,
      magnitude = r$magnitude,
      stringsAsFactors = FALSE
    )
  }))
  rownames(cpf_df) <- cpf_df$cell

  ## Attach metadata colour variable
  if (point_colour %in% colnames(seurat_obj@meta.data)) {
    cpf_df[[point_colour]] <- as.character(
      seurat_obj@meta.data[cpf_df$cell, point_colour])
  }

  ## ── Build ggplot ───────────────────────────────────────────────────────────
  bg_df <- data.frame(
    x = embedding_coords[all_cells, 1],
    y = embedding_coords[all_cells, 2]
  )
  if (point_colour %in% colnames(seurat_obj@meta.data)) {
    bg_df[[point_colour]] <- as.character(
      seurat_obj@meta.data[all_cells, point_colour])
    colour_aes <- ggplot2::aes(colour = .data[[point_colour]])
  } else {
    colour_aes <- ggplot2::aes()
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data    = bg_df,
      mapping = ggplot2::aes(x = x, y = y),
      colour  = "grey80",
      size    = point_size,
      alpha   = 0.6
    ) +
    ggplot2::geom_segment(
      data    = cpf_df,
      mapping = ggplot2::aes(
        x    = x, y    = y,
        xend = x + dx, yend = y + dy
      ),
      arrow   = grid::arrow(length = grid::unit(0.08, "cm"), type = "closed"),
      alpha   = arrow_alpha,
      colour  = "steelblue",
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = title,
      x     = "UMAP 1",
      y     = "UMAP 2"
    ) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text   = ggplot2::element_text(face = "bold"),
      axis.title  = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )

  return(list(cpf_df = cpf_df, plot = p, hvg = hvg))
}


## ── Internal: CPF for one cell ────────────────────────────────────────────────
#' @keywords internal
.cpf_one_cell <- function(cell,
                           gas_hvg,
                           rna_hvg,
                           embedding,
                           knn_k) {
  tryCatch({
    ## Spearman correlation between this cell's scGAS and all cells' RNA
    gas_vec   <- gas_hvg[cell, ]
    cors      <- apply(rna_hvg, 1,
                       function(r) stats::cor(gas_vec, r,
                                              method = "spearman"))
    cors[cell] <- -Inf   # exclude self

    ## Top-k RNA neighbours
    top_k <- names(sort(cors, decreasing = TRUE))[seq_len(knn_k)]

    ## Centroid of those neighbours in embedding space
    target_xy <- colMeans(embedding[top_k, , drop = FALSE])
    self_xy   <- embedding[cell, ]
    delta     <- target_xy - self_xy

    list(
      dx        = delta[1],
      dy        = delta[2],
      magnitude = sqrt(sum(delta^2))
    )
  }, error = function(e) NULL)
}
