#' Write a Cas-OFFinder input file
#'
#' @param protospacer Character vector of protospacers.
#' @param reference Path to the reference FASTA or directory that Cas-OFFinder
#'   will search.
#' @param file Output path.
#' @param pam IUPAC PAM motif.
#' @param mismatches Integer mismatch allowance.
#' @param pam_side `"3prime"` or `"5prime"`.
#' @return `file`, invisibly.
#' @examples
#' f <- tempfile()
#' write_casoffinder_input(c("ACGTACGTACGTACGTACGT"), "panel.fa", f)
#' readLines(f)
#' @export
write_casoffinder_input <- function(protospacer, reference, file, pam = "NGG",
                                    mismatches = 3L,
                                    pam_side = c("3prime", "5prime")) {
  pam_side <- match.arg(pam_side)
  protospacer <- unique(toupper(protospacer))
  if (!length(protospacer)) stop("`protospacer` is empty.", call. = FALSE)
  gl <- unique(nchar(protospacer))
  if (length(gl) != 1L)
    stop("All protospacers must share one length.", call. = FALSE)
  npad <- strrep("N", nchar(pam))
  pattern <- if (identical(pam_side, "3prime"))
    paste0(strrep("N", gl), toupper(pam)) else paste0(toupper(pam), strrep("N", gl))
  query <- if (identical(pam_side, "3prime"))
    paste0(protospacer, npad) else paste0(npad, protospacer)
  writeLines(c(reference, pattern,
               paste(query, as.integer(mismatches))), con = file)
  invisible(file)
}

#' Read Cas-OFFinder output
#'
#' Parses the default tab-separated output. Cas-OFFinder writes mismatched
#' bases in lower case; both the verbatim string and an upper-cased copy are
#' returned, and mismatch positions are resolved at parse time, so downstream
#' joins never depend on letter case.
#'
#' @param file Path to the output file.
#' @param pam_side `"3prime"` or `"5prime"`.
#' @return A data frame with columns `query`, `protospacer`, `allele`,
#'   `position`, `dna`, `dna_upper`, `strand`, `mismatch` and
#'   `mismatch_positions` (a list column).
#' @export
read_casoffinder <- function(file, pam_side = c("3prime", "5prime")) {
  pam_side <- match.arg(pam_side)
  cn <- c("query", "allele", "position", "dna", "strand", "mismatch")
  if (!file.exists(file)) stop("No such file: ", file, call. = FALSE)
  if (!length(readLines(file, n = 1L, warn = FALSE))) {
    out <- data.frame(query = character(), protospacer = character(),
                      allele = character(), position = integer(),
                      dna = character(), dna_upper = character(),
                      strand = character(), mismatch = integer(),
                      stringsAsFactors = FALSE)
    out$mismatch_positions <- list()
    return(out)
  }
  d <- utils::read.delim(file, header = FALSE, col.names = cn,
                         stringsAsFactors = FALSE)
  d$protospacer <- if (identical(pam_side, "3prime"))
    sub("N+$", "", d$query) else sub("^N+", "", d$query)
  d$dna_upper <- toupper(d$dna)
  d$mismatch_positions <- mismatch_positions(d$dna)
  d[, c("query", "protospacer", "allele", "position", "dna", "dna_upper",
        "strand", "mismatch", "mismatch_positions")]
}

#' Positions of lower-cased mismatch bases
#'
#' @param dna Character vector as written by Cas-OFFinder.
#' @return A list of integer vectors.
#' @examples
#' mismatch_positions(c("ACGTacGT", "ACGTACGT"))
#' @export
mismatch_positions <- function(dna) {
  lapply(dna, function(s) {
    if (is.na(s)) return(integer())
    which(strsplit(s, NULL)[[1]] %in% c("a", "c", "g", "t"))
  })
}

#' @noRd
.locus <- function(allele) sub("\\*.*$", "", allele)

#' Off-target counts per protospacer and locus
#'
#' @param hits A data frame with `protospacer` and `allele` columns, as
#'   returned by [read_casoffinder()].
#' @param universe Optional character vector naming every allele in the
#'   searched panel. Supplying it adds the denominator and percentage columns.
#' @return A data frame with one row per protospacer and locus.
#' @examples
#' hits <- data.frame(protospacer = "ACGT", allele = c("A*01:01", "B*07:02"))
#' offtarget_summary(hits, universe = c("A*01:01", "A*02:01", "B*07:02"))
#' @export
offtarget_summary <- function(hits, universe = NULL) {
  if (!all(c("protospacer", "allele") %in% names(hits)))
    stop("`hits` needs `protospacer` and `allele` columns.", call. = FALSE)
  u <- unique(hits[, c("protospacer", "allele")])
  u$locus <- .locus(u$allele)
  agg <- stats::aggregate(allele ~ protospacer + locus, data = u,
                          FUN = function(x) length(unique(x)))
  names(agg)[names(agg) == "allele"] <- "n_alleles"
  if (!is.null(universe)) {
    tab <- table(.locus(unique(universe)))
    agg$n_universe <- as.integer(tab[agg$locus])
    agg$pct_universe <- 100 * agg$n_alleles / agg$n_universe
  }
  agg <- agg[order(agg$protospacer, agg$locus), , drop = FALSE]
  rownames(agg) <- NULL
  agg
}

#' Off-target totals per protospacer
#'
#' @inheritParams offtarget_summary
#' @return A data frame with one row per protospacer.
#' @examples
#' hits <- data.frame(protospacer = "ACGT", allele = c("A*01:01", "B*07:02"))
#' offtarget_totals(hits, universe = c("A*01:01", "A*02:01", "B*07:02"))
#' @export
offtarget_totals <- function(hits, universe = NULL) {
  if (!all(c("protospacer", "allele") %in% names(hits)))
    stop("`hits` needs `protospacer` and `allele` columns.", call. = FALSE)
  u <- unique(hits[, c("protospacer", "allele")])
  out <- stats::aggregate(allele ~ protospacer, data = u,
                          FUN = function(x) length(unique(x)))
  names(out)[names(out) == "allele"] <- "n_alleles"
  if (!is.null(universe)) {
    out$n_universe <- length(unique(universe))
    out$pct_universe <- 100 * out$n_alleles / out$n_universe
  }
  out <- out[order(out$n_alleles), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Drop guides that hit a protected allele
#'
#' @param guides A data frame of guides.
#' @param hits Off-target hits, as returned by [read_casoffinder()].
#' @param pattern Regular expression matched against `hits$allele`, e.g.
#'   `"^B\\*07:02"`.
#' @param protospacer_col Name of the protospacer column in `guides`.
#' @return `guides`, with every protospacer that hits a matching allele removed.
#' @examples
#' guides <- data.frame(protospacer = c("AAA", "CCC"))
#' hits <- data.frame(protospacer = "AAA", allele = "B*07:02")
#' preserving_guides(guides, hits, "^B\\*07:02")
#' @export
preserving_guides <- function(guides, hits, pattern,
                              protospacer_col = "protospacer") {
  if (!protospacer_col %in% names(guides))
    stop("Column '", protospacer_col, "' not found in `guides`.", call. = FALSE)
  bad <- unique(hits$protospacer[grepl(pattern, hits$allele)])
  out <- guides[!guides[[protospacer_col]] %in% bad, , drop = FALSE]
  rownames(out) <- NULL
  out
}
