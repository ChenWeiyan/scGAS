#' Train Per-Gene Lasso Models and Compute Metacell-Level GAS
#'
#' For each gene in \code{Sig.ass.l} that has enough cCRE-Gene associations
#' in the current dataset, this function trains a Lasso regression model using
#' ENCODE bulk DNase-seq and RNA-seq data as training examples.  It then
#' predicts the Gene Activation Score (GAS) at the Metacell level by applying
#' the trained model to the aggregated scATAC-seq signal.
#'
#' Internally the function:
#' \enumerate{
#'   \item Identifies cCREs shared between the reference association list and
#'     the dataset.
#'   \item Aggregates raw counts across cells within each cluster (Metacell).
#'   \item Library-size normalises and log2-transforms the Metacell matrix
#'     using the same scaling factor as the ENCODE training data.
#'   \item Trains a Lasso (\code{alpha = 1}) model per gene with
#'     \code{glmnetUtils::glmnet} and selects lambda at 95% of maximum
#'     deviance explained.
#'   \item Predicts GAS for each Metacell and zero-floors the result.
#' }
#'
#' @param seurat_obj A preprocessed \code{SeuratObject} returned by
#'   \code{\link{scgas_preprocess}}.  Must contain an \code{ATAC} assay and
#'   \code{seurat_clusters} metadata.
#' @param sig_associations Named list of per-gene association objects as
#'   produced by the reference map construction step.  Each element must
#'   contain a \code{Sig.df} data frame with a \code{CRE.sig} column listing
#'   significant cCRE names for that gene.
#' @param encode_rna_lognorm Matrix (genes x samples) of log-normalised
#'   RNA-seq values from ENCODE bulk tissues used for model training.
#' @param encode_dnase_lognorm Matrix (cCREs x samples) of log-normalised
#'   DNase-seq values from the same ENCODE bulk tissues.
#' @param cre_bed_path Path to the reference cCRE BED file, used to match
#'   cCRE names to the \code{SeuratObject} features.
#' @param encode_lib_size Integer.  Total number of cCRE regions in the ENCODE
#'   reference (used to compute the library-size scaling factor).  Defaults to
#'   \code{2081250L}.
#' @param min_associations Integer.  Genes with fewer than this many cCRE
#'   associations in the dataset are skipped.  Passing \code{NULL} (default)
#'   uses the median association count across all genes as the threshold.
#' @param lasso_alpha Numeric elastic-net mixing parameter.  \code{1} (default)
#'   gives pure Lasso; \code{0} gives Ridge.
#' @param deviance_quantile Numeric in (0, 1].  The lambda selected for
#'   prediction is the smallest lambda such that \code{dev.ratio} is at least
#'   \code{deviance_quantile * max(dev.ratio)}.  Default \code{0.95}.
#' @param n_cores Number of cores for \code{parallel::mclapply}.  Default
#'   \code{1}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return A named list with one element per modelled gene, each containing:
#'   \describe{
#'     \item{\code{glmnetfit}}{Fitted \code{glmnet} object.}
#'     \item{\code{lambda_choose}}{Selected regularisation parameter.}
#'     \item{\code{GAS}}{Named numeric vector of Metacell-level GAS values
#'       (zero-floored, same length as number of clusters).}
#'     \item{\code{cCRE_used}}{Character vector of cCRE names used in the
#'       model.}
#'   }
#'   An attribute \code{"ATAC_CL_lognorm"} on the returned list carries the
#'   normalised Metacell ATAC matrix for downstream use by
#'   \code{\link{scgas_compute}}.
#'
#' @examples
#' \dontrun{
#' models <- scgas_train_models(
#'   seurat_obj          = obj,
#'   sig_associations    = sig_ass_list,
#'   encode_rna_lognorm  = rna_mat,
#'   encode_dnase_lognorm = dnase_mat,
#'   cre_bed_path        = "hg19_500bp_CRE.bed",
#'   n_cores             = 8
#' )
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges findOverlaps sort
#' @importFrom Matrix rowSums
#' @importFrom MatrixGenerics rowSums2
#' @importFrom parallel mclapply
#' @importFrom glmnetUtils glmnet
#' @importFrom stats predict median
#' @export
scgas_train_models <- function(seurat_obj,
                               sig_associations,
                               encode_rna_lognorm,
                               encode_dnase_lognorm,
                               cre_bed_path,
                               encode_lib_size     = 2081250L,
                               min_associations    = NULL,
                               lasso_alpha         = 1,
                               deviance_quantile   = 0.95,
                               n_cores             = 1L,
                               verbose             = TRUE) {

  ## ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    is.list(sig_associations), !is.null(names(sig_associations)),
    is.matrix(encode_rna_lognorm) || is.data.frame(encode_rna_lognorm),
    is.matrix(encode_dnase_lognorm),
    file.exists(cre_bed_path),
    is.numeric(encode_lib_size), encode_lib_size > 0,
    is.numeric(lasso_alpha), lasso_alpha >= 0, lasso_alpha <= 1,
    is.numeric(deviance_quantile), deviance_quantile > 0,
    deviance_quantile <= 1
  )

  ## ── Recover cCRE GRanges from the Seurat object ───────────────────────────
  if (verbose) message("[scGAS] Matching cCRE names from reference BED")
  ref_cre_gr <- rtracklayer::import(cre_bed_path)
  obj_ranges <- Signac::GetAssayData(seurat_obj, assay = "ATAC",
                                     slot = "ranges")
  if (is.null(obj_ranges))
    obj_ranges <- seurat_obj@assays$ATAC@ranges

  ov <- GenomicRanges::findOverlaps(ref_cre_gr, obj_ranges)
  obj_ranges$name <- ref_cre_gr[S4Vectors::queryHits(ov)]$name
  names(obj_ranges) <- obj_ranges$name

  ## cCRE names that appear in the association list
  all_ass_cre <- unique(unlist(lapply(sig_associations,
                                      function(x) x$Sig.df$CRE.sig)))
  all_ass_gr  <- GenomicRanges::sort(
    obj_ranges[which(obj_ranges$name %in% all_ass_cre)]
  )
  all_ass_names <- paste(
    as.character(GenomicRanges::seqnames(all_ass_gr)),
    GenomicRanges::start(all_ass_gr),
    GenomicRanges::end(all_ass_gr),
    sep = "-"
  )

  ## ── Metacell aggregation ───────────────────────────────────────────────────
  if (verbose) message("[scGAS] Aggregating counts into Metacells")
  cl_labels <- seurat_obj$seurat_clusters
  cl_levels <- levels(as.factor(cl_labels))

  ## Subset the count matrix to relevant cCREs
  count_mat <- seurat_obj@assays$ATAC@counts
  ## Keep only the features that match our association cCREs
  feat_inter <- intersect(rownames(count_mat), all_ass_names)
  count_mat  <- count_mat[feat_inter, , drop = FALSE]

  ## Build the Metacell aggregate matrix
  atac_cl_mat <- matrix(
    0L,
    nrow = length(feat_inter),
    ncol = length(cl_levels),
    dimnames = list(feat_inter, paste0("Cluster-", cl_levels))
  )
  for (cl in cl_levels) {
    sel <- which(cl_labels == cl)
    atac_cl_mat[, paste0("Cluster-", cl)] <-
      MatrixGenerics::rowSums2(count_mat[, sel, drop = FALSE])
  }

  ## Scaling factor: ratio of cCRE numbers between dataset and reference
  sf <- floor(encode_lib_size *
                length(feat_inter) / nrow(encode_dnase_lognorm))
  if (verbose) message("[scGAS] Library-size scaling factor: ", sf)

  atac_cl_norm   <- t(t(atac_cl_mat) / colSums(atac_cl_mat)) * sf
  atac_cl_lognorm <- log2(atac_cl_norm + 1)

  ## Map cCRE Seurat-style names back to cCRE BED names
  name_map <- setNames(all_ass_gr$name, all_ass_names)
  rownames(atac_cl_lognorm) <- name_map[rownames(atac_cl_lognorm)]

  ## ── Gene filtering ────────────────────────────────────────────────────────
  ass_counts <- vapply(sig_associations, function(x) {
    length(intersect(all_ass_gr$name, x$Sig.df$CRE.sig))
  }, integer(1L))

  thr <- if (is.null(min_associations)) stats::median(ass_counts)
         else as.integer(min_associations)
  sel_genes <- names(which(ass_counts > thr))

  if (length(sel_genes) == 0)
    stop("[scGAS] No genes passed the association threshold (",
         thr, "). Try lowering min_associations.")
  if (verbose) message("[scGAS] Training models for ", length(sel_genes),
                       " genes (threshold = ", thr, " associations)")

  ## ── Parallel model training ───────────────────────────────────────────────
  trained_models <- parallel::mclapply(
    X        = as.list(sel_genes),
    FUN      = .train_one_gene,
    sig_associations   = sig_associations,
    all_ass_gr         = all_ass_gr,
    encode_rna_lognorm = encode_rna_lognorm,
    encode_dnase_lognorm = encode_dnase_lognorm,
    atac_cl_lognorm    = atac_cl_lognorm,
    lasso_alpha        = lasso_alpha,
    deviance_quantile  = deviance_quantile,
    mc.cores = n_cores
  )
  names(trained_models) <- sel_genes

  ## Remove genes where training failed (NULL entries)
  failed <- vapply(trained_models, is.null, logical(1L))
  if (any(failed)) {
    warning("[scGAS] Training failed for ", sum(failed), " gene(s): ",
            paste(sel_genes[failed], collapse = ", "))
    trained_models <- trained_models[!failed]
  }

  ## Attach the Metacell matrix as an attribute for downstream functions
  attr(trained_models, "ATAC_CL_lognorm") <- atac_cl_lognorm
  attr(trained_models, "cl_levels")       <- cl_levels
  attr(trained_models, "all_ass_gr")      <- all_ass_gr

  if (verbose) message("[scGAS] Model training complete for ",
                       length(trained_models), " genes.")
  return(trained_models)
}


## ── Internal: train one gene ──────────────────────────────────────────────
#' @keywords internal
.train_one_gene <- function(gene_name,
                            sig_associations,
                            all_ass_gr,
                            encode_rna_lognorm,
                            encode_dnase_lognorm,
                            atac_cl_lognorm,
                            lasso_alpha,
                            deviance_quantile) {
  tryCatch({
    ## cCREs for this gene that are present in the dataset
    gene_cres <- sig_associations[[gene_name]]$Sig.df$CRE.sig
    cre_v     <- intersect(all_ass_gr$name, gene_cres)
    if (length(cre_v) == 0) return(NULL)

    ## Build training data frame
    train_df <- data.frame(
      RNA = encode_rna_lognorm[gene_name, ],
      t(encode_dnase_lognorm[cre_v, , drop = FALSE])
    )

    ## Lasso fit
    fit <- glmnetUtils::glmnet(
      RNA ~ .,
      data   = train_df,
      family = "gaussian",
      alpha  = lasso_alpha
    )
    class(fit) <- "glmnet"   # ensure predict.glmnet dispatch

    ## Lambda selection: first index where deviance >= quantile * max deviance
    which_max      <- which(
      fit$dev.ratio >= max(fit$dev.ratio) * deviance_quantile
    )[1]
    lambda_choose  <- fit$lambda[which_max]

    ## Metacell-level prediction
    newx <- t(atac_cl_lognorm[cre_v, , drop = FALSE])
    gas  <- as.vector(
      stats::predict(fit, newx = newx, type = "response", s = lambda_choose)
    )
    gas  <- gas - min(gas)   # zero-floor

    list(
      glmnetfit     = fit,
      lambda_choose = lambda_choose,
      GAS           = gas,
      cCRE_used     = cre_v
    )
  }, error = function(e) {
    NULL
  })
}
