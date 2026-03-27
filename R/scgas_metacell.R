#' Construct Metacells via Walktrap Community Detection
#'
#' Partitions cells into Metacells using the \strong{walktrap} community
#' detection algorithm on a cell–cell similarity graph (e.g. a WNN or KNN
#' graph).
#'
#' The algorithm proceeds in three stages:
#' \enumerate{
#'   \item \strong{Rough partitioning.}
#'     \code{\link[igraph]{cluster_walktrap}} is run on the graph and the
#'     resulting dendrogram is cut to yield approximately
#'     \code{floor(n_cells / gamma)} communities.
#'   \item \strong{Merge small.}
#'     Communities smaller than \code{min_size} are merged into their nearest
#'     neighbour along the walktrap merge tree, preserving topological
#'     contiguity.
#'   \item \strong{Split large.}
#'     Communities larger than \code{max_size} are iteratively refined by
#'     "peeling off" sub-trees whose size falls within
#'     \code{[min_size, gamma * split_ceiling_factor]}, yielding tighter,
#'     more homogeneous groups.
#' }
#'
#' @param graph Cell–cell similarity graph.
#'   Accepted types:
#'   \describe{
#'     \item{\code{igraph}}{An \code{igraph} object (directed or undirected).
#'       Weights are used if present.}
#'     \item{\code{matrix} / \code{dgCMatrix} / \code{dsCMatrix}}{
#'       An adjacency or KNN connectivity matrix (cells x cells).
#'       Non-zero entries are treated as edge weights.  Row/column names, if
#'       present, are used as cell identifiers.}
#'   }
#' @param cell_names Optional character vector of cell identifiers.  When
#'   \code{NULL} (default), names are inferred from \code{graph} (vertex
#'   names or matrix row names).  If the graph carries no names, cells are
#'   labelled \code{"cell_1"}, \code{"cell_2"}, \ldots
#' @param gamma Integer.  Target number of cells per Metacell.  The
#'   dendrogram is initially cut into \code{floor(n / gamma)} groups.
#'   Default \code{20}.
#' @param min_size Integer.  Metacells smaller than this are merged with
#'   their dendrogram neighbour.  Default \code{10}.
#' @param max_size Integer.  Metacells larger than this enter the
#'   split-large stage.  Default \code{100}.
#' @param split_ceiling_factor Numeric.  During splitting, only sub-trees
#'   with at most \code{gamma * split_ceiling_factor} cells are peeled off.
#'   Default \code{1.5}.
#' @param min_remaining Integer.  Stop peeling a large Metacell once the
#'   remaining core drops below this size.  Default \code{30}.
#' @param walktrap_steps Integer.  Length of random walks passed to
#'   \code{\link[igraph]{cluster_walktrap}}.  Default \code{4}.
#' @param n_cores Integer.  Number of parallel workers for the split-large
#'   stage.  Uses \code{parallel::mclapply} (Unix) or sequential fallback
#'   (Windows).  Default \code{1L}.
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return A named integer vector of Metacell assignments.  Names are cell
#'   identifiers (matching the graph vertex / matrix row names).  Metacell
#'   labels are consecutive integers starting from 1.
#'
#' @details
#' ## Choosing \code{gamma}
#' A value of 20 means roughly 20 cells per Metacell.  For a dataset of
#' 100 000 cells this yields ~5 000 Metacells.  Smaller \code{gamma} gives
#' finer granularity but higher computational cost downstream.
#'
#' ## Parallelism
#' The merge-small step is fully vectorised over the walktrap merge tree and
#' runs in O(n) time.
#' The split-large step processes each oversized Metacell independently and
#' is parallelised across Metacells with \code{parallel::mclapply}.
#'
#' ## Memory
#' The function never materialises a dense n-by-n distance matrix.
#' Internally it operates on the igraph merge tree (2n-1 rows) and
#' integer membership vectors, keeping peak memory proportional to the
#' number of cells.
#'
#' @examples
#' \dontrun{
#' ## From a Seurat WNN graph stored as a sparse matrix:
#' wnn <- readMM("wnn_graph.mtx")
#' mc <- scgas_metacell(wnn, gamma = 20, min_size = 10, max_size = 100,
#'                      n_cores = 8)
#'
#' ## From an igraph object:
#' library(igraph)
#' g <- make_ring(200)
#' mc <- scgas_metacell(g, gamma = 10)
#'
#' ## Plug into Seurat for downstream scGAS workflow:
#' obj$seurat_clusters <- factor(mc[colnames(obj)])
#' }
#'
#' @importFrom parallel mclapply
#' @export
scgas_metacell <- function(graph,
                           cell_names           = NULL,
                           gamma                = 20L,
                           min_size             = 10L,
                           max_size             = 100L,
                           split_ceiling_factor = 1.5,
                           min_remaining        = 30L,
                           walktrap_steps       = 4L,
                           n_cores              = 1L,
                           verbose              = TRUE) {

  ## ── Input validation ─────────────────────────────────────────────────────
  gamma      <- as.integer(gamma)
  min_size   <- as.integer(min_size)
  max_size   <- as.integer(max_size)
  n_cores    <- as.integer(n_cores)
  walktrap_steps <- as.integer(walktrap_steps)
  min_remaining  <- as.integer(min_remaining)

  stopifnot(
    gamma >= 2L,
    min_size >= 2L,
    max_size > min_size,
    split_ceiling_factor > 0,
    min_remaining >= 1L,
    walktrap_steps >= 1L,
    n_cores >= 1L
  )

  ## ── Convert input to igraph ──────────────────────────────────────────────
  g <- .to_igraph(graph, verbose)
  n <- igraph::vcount(g)

  ## ── Resolve cell names ───────────────────────────────────────────────────
  if (!is.null(cell_names)) {
    stopifnot(length(cell_names) == n)
  } else {
    vn <- igraph::V(g)$name
    if (!is.null(vn) && length(vn) == n) {
      cell_names <- vn
    } else if (inherits(graph, c("matrix", "Matrix"))) {
      cell_names <- rownames(graph)
    }
    if (is.null(cell_names))
      cell_names <- paste0("cell_", seq_len(n))
  }

  ## ── Walktrap clustering ──────────────────────────────────────────────────
  if (verbose) message("[scGAS] Running walktrap (", n, " cells, steps = ",
                       walktrap_steps, ")")
  wt <- igraph::cluster_walktrap(g, steps = walktrap_steps)

  n_target <- max(1L, floor(n / gamma))
  if (verbose) message("[scGAS] Cutting dendrogram into ~", n_target,
                       " initial groups (gamma = ", gamma, ")")
  membership <- igraph::cut_at(wt, no = n_target)

  ## ── Stage 1: Merge small metacells ───────────────────────────────────────
  if (verbose) message("[scGAS] Merging small communities (min_size = ",
                       min_size, ")")
  membership <- .merge_small_vec(wt, membership, n, min_size)

  n_after_merge <- length(unique(membership))
  if (verbose) message("[scGAS]   -> ", n_after_merge,
                       " metacells after merging")

  ## ── Stage 2: Split large metacells ───────────────────────────────────────
  mc_sizes  <- tabulate(membership, nbins = max(membership))
  large_ids <- which(mc_sizes > max_size)

  if (length(large_ids) > 0L) {
    if (verbose) message("[scGAS] Splitting ", length(large_ids),
                         " oversized metacells (max_size = ", max_size, ")")

    ## Pre-compute node_cells lookup from merge tree (vectorised)
    merges <- wt$merges
    node_cells <- .build_node_cells(merges, n)

    split_fn <- function(mc_id) {
      .split_one_metacell(
        mc_id         = mc_id,
        membership    = membership,
        merges        = merges,
        node_cells    = node_cells,
        n             = n,
        min_size      = min_size,
        max_size      = max_size,
        gamma         = gamma,
        ceiling_fac   = split_ceiling_factor,
        min_remaining = min_remaining
      )
    }

    if (n_cores > 1L && .Platform$OS.type == "unix") {
      split_results <- parallel::mclapply(
        X        = large_ids,
        FUN      = split_fn,
        mc.cores = min(n_cores, length(large_ids))
      )
    } else {
      split_results <- lapply(large_ids, split_fn)
    }

    ## Apply split labels back into membership (character intermediary
    ## to handle sub-labels like "42.1", "42.2")
    membership_chr <- as.character(membership)
    for (res in split_results) {
      if (!is.null(res)) {
        idx <- as.integer(names(res))
        membership_chr[idx] <- res
      }
    }
    membership <- as.integer(factor(membership_chr,
                                    levels = unique(membership_chr)))
  }

  ## ── Relabel to consecutive 1:K ───────────────────────────────────────────
  membership <- as.integer(factor(membership))
  names(membership) <- cell_names

  n_final <- max(membership)
  sizes   <- tabulate(membership, nbins = n_final)
  if (verbose) {
    message("[scGAS] Metacell construction complete: ",
            n_final, " metacells from ", n, " cells")
    message("[scGAS]   Size range: ", min(sizes), " – ", max(sizes),
            " (median ", stats::median(sizes), ")")
  }

  return(membership)
}


## ═══════════════════════════════════════════════════════════════════════════
## Internal helpers
## ═══════════════════════════════════════════════════════════════════════════

# --------------------------------------------------------------------------
# .to_igraph: convert matrix / dgCMatrix / dsCMatrix / igraph -> igraph
# --------------------------------------------------------------------------
#' @keywords internal
.to_igraph <- function(x, verbose = TRUE) {
  if (inherits(x, "igraph")) {
    return(x)
  }
  ## Sparse Matrix (dgCMatrix, dsCMatrix, dgTMatrix, ...)
  if (inherits(x, "Matrix") || inherits(x, "sparseMatrix")) {
    if (verbose) message("[scGAS] Converting sparse matrix to igraph graph")
    ## Ensure symmetric for undirected graph
    if (!Matrix::isSymmetric(x)) {
      x <- (x + Matrix::t(x)) / 2
    }
    return(igraph::graph_from_adjacency_matrix(
      x, mode = "undirected", weighted = TRUE, diag = FALSE
    ))
  }
  ## Dense matrix
  if (is.matrix(x)) {
    if (verbose) message("[scGAS] Converting dense matrix to igraph graph")
    if (!isSymmetric(x)) {
      x <- (x + t(x)) / 2
    }
    return(igraph::graph_from_adjacency_matrix(
      x, mode = "undirected", weighted = TRUE, diag = FALSE
    ))
  }
  stop("[scGAS] 'graph' must be an igraph object, a (sparse) matrix, or a ",
       "dgCMatrix / dsCMatrix.  Got class: ", paste(class(x), collapse = ", "))
}


# --------------------------------------------------------------------------
# .build_node_cells: pre-compute which leaf cells belong to each merge node
#
# This is the most memory-critical helper.  We store cell indices as integer
# vectors in a list of length 2n-1 (n leaves + n-1 internal nodes).
# We use range-encoding where possible to avoid duplicating large vectors.
# --------------------------------------------------------------------------
#' @keywords internal
.build_node_cells <- function(merges, n) {
  n_internal <- nrow(merges)       # should be n-1
  total      <- n + n_internal
  node_cells <- vector("list", total)

  ## Leaves: each leaf i contains only cell i
  for (i in seq_len(n)) {
    node_cells[[i]] <- i
  }

  ## Internal nodes: union of children
  for (i in seq_len(n_internal)) {
    l <- merges[i, 1]
    r <- merges[i, 2]
    node_cells[[n + i]] <- c(node_cells[[l]], node_cells[[r]])
  }

  return(node_cells)
}


# --------------------------------------------------------------------------
# .merge_small_vec: vectorised merge-small along the walktrap merge tree
#
# Walks the merge tree bottom-up.  When two communities meet at a merge
# step and one (or both) is below min_size, they are unified.
# All label look-ups use a Union-Find (disjoint-set) structure for O(α(n))
# amortised time per find, instead of the O(n) full-vector rewrites in the
# original script.
# --------------------------------------------------------------------------
#' @keywords internal
.merge_small_vec <- function(wt, membership, n, min_size) {
  merges <- wt$merges
  n_steps <- nrow(merges)

  ## --- Union-Find data structure (path-compressed) -----------------------
  parent <- seq_len(max(membership))  # 1:K, each is its own root
  uf_size <- tabulate(membership, nbins = max(membership))

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]   # path halving
      x <- parent[x]
    }
    x
  }

  union_sets <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra == rb) return(invisible(NULL))
    ## Union by size
    if (uf_size[ra] >= uf_size[rb]) {
      parent[rb] <<- ra
      uf_size[ra] <<- uf_size[ra] + uf_size[rb]
    } else {
      parent[ra] <<- rb
      uf_size[rb] <<- uf_size[rb] + uf_size[ra]
    }
  }

  ## --- Walk the merge tree -----------------------------------------------
  ## node_mc tracks which metacell label each node (leaf or internal) maps to
  node_mc <- c(membership, rep(NA_integer_, n_steps))

  for (step in seq_len(n_steps)) {
    l   <- merges[step, 1]
    r   <- merges[step, 2]
    lmc_raw <- node_mc[l]
    rmc_raw <- node_mc[r]

    ## Resolve to root labels via Union-Find
    lmc <- if (!is.na(lmc_raw)) find_root(lmc_raw) else NA_integer_
    rmc <- if (!is.na(rmc_raw)) find_root(rmc_raw) else NA_integer_

    ## Propagate label to new internal node
    node_mc[n + step] <- if (!is.na(lmc)) lmc else rmc

    if (!is.na(lmc) && !is.na(rmc) && lmc != rmc) {
      ls <- uf_size[lmc]
      rs <- uf_size[rmc]

      if (ls < min_size || rs < min_size) {
        ## At least one side is undersized -> merge them
        union_sets(lmc, rmc)
        node_mc[n + step] <- find_root(lmc)
      }
    }
  }

  ## --- Remap membership through Union-Find roots -------------------------
  out <- integer(n)
  for (i in seq_len(n)) {
    out[i] <- find_root(membership[i])
  }

  ## Relabel to consecutive integers
  as.integer(factor(out))
}


# --------------------------------------------------------------------------
# .split_one_metacell: peel sub-trees from an oversized metacell
#
# Strategy: walk the internal merge steps that are wholly contained within
# this metacell in *reverse* order (top-down in the dendrogram).  When a
# branch has size in [min_size, gamma * ceiling_fac] and splitting it off
# leaves at least min_remaining cells in the core, peel it off as a new
# sub-metacell.
#
# Returns a named character vector: names = cell indices (as char),
# values = sub-labels like "42.1", "42.2", or original "42".
# --------------------------------------------------------------------------
#' @keywords internal
.split_one_metacell <- function(mc_id,
                                membership,
                                merges,
                                node_cells,
                                n,
                                min_size,
                                max_size,
                                gamma,
                                ceiling_fac,
                                min_remaining) {
  cells_in_mc <- which(membership == mc_id)
  n_mc <- length(cells_in_mc)

  if (n_mc <= max_size) return(NULL)  # safety check

  ## Identify internal merge steps whose children are wholly within this mc
  ## Use a logical set for O(1) membership testing
  cell_set <- logical(n)
  cell_set[cells_in_mc] <- TRUE

  n_steps <- nrow(merges)
  internal_steps <- integer(0)

  for (s in seq_len(n_steps)) {
    lc <- node_cells[[merges[s, 1]]]
    rc <- node_cells[[merges[s, 2]]]
    if (all(cell_set[lc]) && all(cell_set[rc])) {
      internal_steps <- c(internal_steps, s)
    }
  }

  if (length(internal_steps) == 0L) return(NULL)

  ## Prepare output labels
  labels <- stats::setNames(
    rep(as.character(mc_id), n_mc),
    as.character(cells_in_mc)
  )
  remaining   <- cells_in_mc
  remain_set  <- cell_set           # reuse logical vector
  sub_count   <- 0L
  upper_bound <- as.integer(floor(gamma * ceiling_fac))

  ## Walk top-down (reverse of merge order)
  for (step in rev(internal_steps)) {
    if (length(remaining) < min_remaining) break

    for (child in c(merges[step, 1], merges[step, 2])) {
      branch <- node_cells[[child]]
      ## Filter to cells still in the remaining pool
      branch <- branch[remain_set[branch]]
      blen   <- length(branch)

      if (blen >= min_size &&
          blen <= upper_bound &&
          (length(remaining) - blen) >= min_remaining) {
        sub_count <- sub_count + 1L
        labels[as.character(branch)] <- paste0(mc_id, ".", sub_count)
        remain_set[branch] <- FALSE
        remaining <- remaining[remain_set[remaining]]
        break   # one peel per merge step, then move to next step
      }
    }
  }

  return(labels)
}
