#' Preprocess scATAC-seq Data Against Reference cCREs
#'
#' Builds a \code{SeuratObject} with an \code{ATAC} assay whose features are
#' the reference cCREs supplied by the user, runs standard quality control
#' metrics (nucleosome signal, TSS enrichment), and returns a dimensionality-
#' reduced object ready for Metacell construction.
#'
#' The function supports **single-replicate** and **multi-replicate** datasets.
#' When more than one fragment file is provided the per-replicate count
#' matrices are column-bound before all downstream steps.
#'
#' @param fragment_paths Character vector of paths to fragment files
#'   (\code{.tsv.gz}).  One element per replicate.
#' @param cre_bed_path Path to the reference cCRE BED file
#'   (\code{hg19_500bp_CRE.bed} or equivalent).  The file must contain a
#'   \code{name} column that uniquely identifies every 500 bp bin.
#' @param rna_matrix Optional pre-loaded sparse RNA count matrix
#'   (genes x cells).  When supplied the function adds an \code{RNA} assay to
#'   the returned object and runs basic RNA processing (variable features,
#'   scaling, PCA).
#' @param meta_data Optional data frame of per-cell metadata.  Row names must
#'   match cell barcodes.
#' @param annotation \code{GRanges} object of gene annotations compatible with
#'   the genome build of the fragment files.  Obtain with
#'   \code{Signac::GetGRangesFromEnsDb()} and then set
#'   \code{seqlevelsStyle(annotation) <- "UCSC"}.
#' @param cell_barcodes Optional character vector of cell barcodes to retain.
#'   Useful when RNA and ATAC barcodes differ by a replicate prefix.  If
#'   \code{NULL} (default) all barcodes present in the feature matrix are
#'   kept.
#' @param replicate_ids Optional character vector of replicate identifiers
#'   (same length as \code{fragment_paths}).  Used as \code{add.cell.ids} when
#'   merging multi-replicate objects.  Defaults to
#'   \code{"rep1", "rep2", ...}.
#' @param min_cells_per_feature Minimum number of cells a cCRE must be
#'   detected in to be retained.  Default \code{1}.
#' @param min_features_per_cell Minimum number of cCREs detected per cell.
#'   Default \code{1}.
#' @param lsi_dims Integer vector of LSI components to use for UMAP and
#'   neighbour graph construction.  Component 1 is almost always correlated
#'   with sequencing depth and should be excluded.  Default \code{2:30}.
#' @param n_cores Number of cores passed to \code{parallel::mclapply}.
#'   Default \code{1}.
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return A \code{SeuratObject} with:
#'   \describe{
#'     \item{\code{ATAC}}{Chromatin assay containing raw cCRE counts and LSI
#'       reduction.}
#'     \item{\code{RNA}}{(Optional) RNA assay with log-normalised data and PCA
#'       reduction.}
#'     \item{\code{umap}}{3-component UMAP embedding for visualisation and
#'       trajectory analysis.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v75)
#' annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v75)
#' seqlevelsStyle(annotation) <- "UCSC"
#' genome(annotation) <- "hg19"
#'
#' obj <- scgas_preprocess(
#'   fragment_paths = c("rep1_fragments.tsv.gz", "rep2_fragments.tsv.gz"),
#'   cre_bed_path   = "hg19_500bp_CRE.bed",
#'   annotation     = annotation,
#'   replicate_ids  = c("rep1", "rep2"),
#'   lsi_dims       = 2:30,
#'   cluster_resolution = 1.4
#' )
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom Signac CreateFragmentObject FeatureMatrix CreateChromatinAssay
#'   NucleosomeSignal TSSEnrichment FindTopFeatures RunTFIDF RunSVD
#' @importFrom Seurat CreateSeuratObject CreateAssayObject FindNeighbors FindClusters RunUMAP AddMetaData
#' @importFrom Matrix Matrix
#' @importFrom GenomicRanges findOverlaps seqnames
#' @importFrom GenomeInfoDb seqlevels 
#' @export
scgas_preprocess <- function(fragment_paths,
                             cre_bed_path,
                             rna_matrix        = NULL,
                             meta_data         = NULL,
                             annotation        = NULL,
                             cell_barcodes     = NULL,
                             replicate_ids     = NULL,
                             min_cells_per_feature  = 1L,
                             min_features_per_cell  = 1L,
                             lsi_dims          = 2:30,
                             n_cores           = 1L,
                             verbose           = TRUE) {
  
  ## ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    is.character(fragment_paths), length(fragment_paths) >= 1,
    is.character(cre_bed_path),   length(cre_bed_path) == 1,
    file.exists(cre_bed_path)
  )
  for (fp in fragment_paths) {
    if (!file.exists(fp))
      stop("Fragment file not found: ", fp)
  }
  
  n_reps <- length(fragment_paths)
  if (is.null(replicate_ids))
    replicate_ids <- paste0("rep", seq_len(n_reps))
  stopifnot(length(replicate_ids) == n_reps)
  
  ## ── Load reference cCRE ranges ─────────────────────────────────────────────
  if (verbose) message("[scGAS] Loading reference cCRE BED: ", cre_bed_path)
  ref_cre_gr <- rtracklayer::import(cre_bed_path)
  names(ref_cre_gr) <- ref_cre_gr$name
  
  ## ── Harmonise chromosome naming between BED and fragment files ─────────────
  ## Detect whether the first fragment file uses "chr"-prefixed chromosomes by
  ## reading the first non-header line.  The BED is then adjusted to match, so
  ## FeatureMatrix() finds overlapping seqlevels regardless of which convention
  ## the user's genome / fragment files follow.
  ref_cre_gr <- .harmonise_seqlevels(ref_cre_gr, fragment_paths[1], verbose)
  
  ## ── Build per-replicate count matrices ─────────────────────────────────────
  count_list <- vector("list", n_reps)
  frag_list  <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    if (verbose) message("[scGAS] Processing replicate ", replicate_ids[i])
    frag_list[[i]] <- Signac::CreateFragmentObject(fragment_paths[i])
    cm <- Signac::FeatureMatrix(
      fragments = frag_list[[i]],
      features  = ref_cre_gr
    )
    ## Optionally restrict to known barcodes
    if (!is.null(cell_barcodes)) {
      bc_inter <- intersect(colnames(cm), cell_barcodes)
      cm <- cm[, bc_inter, drop = FALSE]
    }
    count_list[[i]] <- cm
  }
  
  ## ── Merge replicates ───────────────────────────────────────────────────────
  if (n_reps == 1) {
    merged_counts <- count_list[[1]]
  } else {
    merged_counts <- do.call(cbind, count_list)
  }
  
  ## Remove all-zero cCREs
  zero_rows <- which(Matrix::rowSums(merged_counts) == 0)
  if (length(zero_rows) > 0) {
    if (verbose) message("[scGAS] Removing ", length(zero_rows),
                         " all-zero cCREs")
    merged_counts <- merged_counts[-zero_rows, , drop = FALSE]
    ref_cre_gr    <- ref_cre_gr[-zero_rows]
  }
  
  ## ── Build SeuratObject per replicate then merge ────────────────────────────
  seurat_list <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    rep_counts <- count_list[[i]]
    ## Align features after zero-row removal
    rep_counts <- rep_counts[rownames(merged_counts), , drop = FALSE]
    
    chrom_assay <- Signac::CreateChromatinAssay(
      fragments    = frag_list[[i]],
      counts       = rep_counts,
      ranges       = ref_cre_gr,
      min.cells    = min_cells_per_feature,
      min.features = min_features_per_cell,
      annotation   = annotation
    )
    seurat_list[[i]] <- Seurat::CreateSeuratObject(
      counts = chrom_assay,
      assay  = "ATAC"
    )
  }
  
  if (n_reps == 1) {
    obj <- seurat_list[[1]]
  } else {
    ## merge() is the correct S3 generic for SeuratObjects in both Seurat v4
    ## and v5.  Passing y as a list and add.cell.ids as a vector of length
    ## n_reps prepends each replicate prefix exactly once across all barcodes.
    obj <- merge(
      x            = seurat_list[[1]],
      y            = seurat_list[-1],   # list of remaining SeuratObjects
      add.cell.ids = replicate_ids      # one prefix per replicate
    )
  }
  
  ## ── Add RNA assay (optional) ───────────────────────────────────────────────
  if (!is.null(rna_matrix)) {
    if (verbose) message("[scGAS] Adding RNA assay")
    rna_matrix <- Matrix::Matrix(rna_matrix, sparse = TRUE)
    ## Align cell names
    rna_matrix <- rna_matrix[, intersect(colnames(rna_matrix),
                                         colnames(obj)), drop = FALSE]
    obj[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_matrix)
    Seurat::DefaultAssay(obj) <- "RNA"
    obj <- Seurat::FindVariableFeatures(obj, verbose = verbose)
    obj <- Seurat::ScaleData(obj, verbose = verbose)
    obj <- Seurat::RunPCA(obj, verbose = verbose)
  }
  
  ## ── Add metadata (optional) ────────────────────────────────────────────────
  if (!is.null(meta_data)) {
    shared_cells <- intersect(rownames(meta_data), colnames(obj))
    obj <- Seurat::AddMetaData(obj, metadata = meta_data[shared_cells, ,
                                                         drop = FALSE])
  }
  
  ## ── QC metrics ────────────────────────────────────────────────────────────
  if (verbose) message("[scGAS] Computing QC metrics")
  Seurat::DefaultAssay(obj) <- "ATAC"
  obj <- Signac::NucleosomeSignal(obj)
  obj <- Signac::TSSEnrichment(obj)
  
  ## ── ATAC dimensionality reduction ─────────────────────────────────────────
  if (verbose) message("[scGAS] Running LSI")
  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q10", verbose = verbose)
  obj <- Signac::RunTFIDF(obj, verbose = verbose)
  obj <- Signac::RunSVD(obj, verbose = verbose)
  
  ## ── 3-component UMAP for trajectory analysis ──────────────────────────────
  if (verbose) message("[scGAS] Running UMAP (3 components)")
  obj <- Seurat::RunUMAP(
    object      = obj,
    reduction   = "lsi",
    dims        = lsi_dims,
    n.components = 3L,
    verbose     = verbose
  )
  
  if (verbose) message("[scGAS] Preprocessing complete. ",
                       ncol(obj), " cells, ",
                       nrow(obj), " cCREs retained.")
  return(obj)
}


## ── Internal: harmonise seqlevels between BED and fragment file ───────────────
##
## Strategy:
##   1. Peek at the first data line of the fragment file (skip comment/header
##      lines starting with "#").
##   2. Determine whether the fragment file uses "chr"-prefixed chromosomes.
##   3. Check what the BED currently uses.
##   4. Strip or add "chr" on the BED GRanges to match the fragment convention.
##
## This means the BED can be delivered in either UCSC ("chr1") or NCBI ("1")
## style and will be silently corrected to match the fragments.
##
#' @keywords internal
.harmonise_seqlevels <- function(gr, fragment_path, verbose = TRUE) {
  
  ## --- Detect fragment chromosome style ------------------------------------
  frag_chr_style <- tryCatch({
    con <- if (grepl("\\.gz$", fragment_path, ignore.case = TRUE))
      gzcon(file(fragment_path, "rb"))
    else
      file(fragment_path, "r")
    on.exit(close(con), add = TRUE)
    
    first_chrom <- NULL
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0) break        # EOF
      if (!startsWith(line, "#")) {       # skip comment / header lines
        first_chrom <- strsplit(line, "\t")[[1]][1]
        break
      }
    }
    if (is.null(first_chrom)) "unknown"
    else if (startsWith(first_chrom, "chr")) "UCSC"
    else "NCBI"
  }, error = function(e) "unknown")
  
  if (frag_chr_style == "unknown") {
    if (verbose)
      message("[scGAS] Could not detect chromosome style from fragment file; ",
              "leaving BED seqlevels unchanged.")
    return(gr)
  }
  
  ## --- Detect BED chromosome style ----------------------------------------
  bed_has_chr <- any(startsWith(
    as.character(GenomicRanges::seqnames(gr)), "chr"
  ))
  bed_style <- if (bed_has_chr) "UCSC" else "NCBI"
  
  if (bed_style == frag_chr_style) {
    if (verbose)
      message("[scGAS] Chromosome style consistent (", frag_chr_style,
              "); no seqlevels adjustment needed.")
    return(gr)
  }
  
  ## --- Adjust BED to match fragment style ----------------------------------
  if (frag_chr_style == "NCBI" && bed_style == "UCSC") {
    ## BED has "chr", fragments do not -> strip "chr" from BED seqlevels
    if (verbose)
      message("[scGAS] Fragment file uses NCBI-style chromosomes (no 'chr'); ",
              "stripping 'chr' prefix from BED seqlevels.")
    GenomicRanges::seqlevels(gr) <-
      gsub("^chr", "", GenomicRanges::seqlevels(gr))
    
  } else {
    ## BED has no "chr", fragments do -> add "chr" to BED seqlevels
    if (verbose)
      message("[scGAS] Fragment file uses UCSC-style chromosomes ('chr'); ",
              "adding 'chr' prefix to BED seqlevels.")
    GenomicRanges::seqlevels(gr) <-
      paste0("chr", GenomicRanges::seqlevels(gr))
  }
  
  return(gr)
}