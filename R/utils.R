#' IUPAC nucleotide codes
#' @noRd
.iupac <- c(A = "A", C = "C", G = "G", T = "T",
            R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC",
            B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT")

#' Reverse complement a DNA sequence
#'
#' @param x Character vector of DNA sequences. IUPAC ambiguity codes are
#'   complemented; case is preserved.
#' @return A character vector of the same length as `x`.
#' @examples
#' revcomp(c("ACGTN", "TTTV"))
#' @export
revcomp <- function(x) {
  if (!is.character(x)) stop("`x` must be a character vector.", call. = FALSE)
  comp <- chartr("ACGTacgtRYKMBVDHrykmbvdhSWNswn",
                 "TGCAtgcaYRMKVBHDyrmkvbhdSWNswn", x)
  vapply(comp, function(s) {
    if (is.na(s)) return(NA_character_)
    paste(rev(strsplit(s, NULL)[[1]]), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

#' Translate an IUPAC motif into a regular expression
#'
#' @param motif A character scalar such as `"NGG"` or `"TTTV"`.
#' @return A character scalar holding a regular expression.
#' @examples
#' iupac_regex("NGG")
#' iupac_regex("TTTV")
#' @export
iupac_regex <- function(motif) {
  if (!is.character(motif) || length(motif) != 1L || is.na(motif) || !nzchar(motif))
    stop("`motif` must be a non-empty character scalar.", call. = FALSE)
  ch <- strsplit(toupper(motif), NULL)[[1]]
  bad <- setdiff(ch, names(.iupac))
  if (length(bad))
    stop("Not an IUPAC code: ", paste(bad, collapse = ", "), call. = FALSE)
  paste0(vapply(ch, function(c1) {
    s <- .iupac[[c1]]
    if (nchar(s) == 1L) s else paste0("[", s, "]")
  }, character(1)), collapse = "")
}

#' Enumerate protospacers adjacent to a PAM
#'
#' Scans both strands of `sequence` for the PAM motif and returns every
#' protospacer of the requested length that sits adjacent to one. Overlapping
#' PAM occurrences are all reported. Coordinates are always given on the plus
#' strand of `sequence` with `start <= end`; the `strand` column records which
#' strand the protospacer is read from.
#'
#' @param sequence Character scalar; the DNA sequence to scan.
#' @param pam IUPAC PAM motif, e.g. `"NGG"` or `"TTTV"`.
#' @param guide_length Integer protospacer length.
#' @param pam_side `"3prime"` for Cas9-like PAMs that follow the protospacer,
#'   `"5prime"` for Cas12a-like PAMs that precede it.
#' @param strands Which strands to scan; defaults to both.
#' @return A data frame with columns `strand`, `protospacer`, `pam`, `start`,
#'   `end` and `pam_start`. Zero rows if nothing matches.
#' @examples
#' find_protospacers("ACGTACGTACGTACGTACGTAGGTT", pam = "NGG", guide_length = 20)
#' @export
find_protospacers <- function(sequence, pam = "NGG", guide_length = 20L,
                              pam_side = c("3prime", "5prime"),
                              strands = c("+", "-")) {
  pam_side <- match.arg(pam_side)
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence))
    stop("`sequence` must be a character scalar.", call. = FALSE)
  guide_length <- as.integer(guide_length)
  if (is.na(guide_length) || guide_length < 1L)
    stop("`guide_length` must be a positive integer.", call. = FALSE)
  strands <- match.arg(strands, c("+", "-"), several.ok = TRUE)

  sequence <- toupper(sequence)
  L  <- nchar(sequence)
  rx <- iupac_regex(pam)
  np <- nchar(pam)
  out <- list()

  for (st in strands) {
    s <- if (identical(st, "+")) sequence else revcomp(sequence)
    hit <- as.integer(gregexpr(paste0("(?=", rx, ")"), s, perl = TRUE)[[1]])
    hit <- hit[hit > 0L]
    if (!length(hit)) next
    if (identical(pam_side, "3prime")) {
      ps <- hit - guide_length; pe <- hit - 1L
    } else {
      ps <- hit + np; pe <- hit + np + guide_length - 1L
    }
    keep <- ps >= 1L & pe <= L & (hit + np - 1L) <= L
    hit <- hit[keep]; ps <- ps[keep]; pe <- pe[keep]
    if (!length(hit)) next
    if (identical(st, "+")) {
      gs <- ps; ge <- pe; pstart <- hit
    } else {
      gs <- L - pe + 1L; ge <- L - ps + 1L; pstart <- L - (hit + np - 1L) + 1L
    }
    out[[st]] <- data.frame(
      strand = st, protospacer = substring(s, ps, pe),
      pam = substring(s, hit, hit + np - 1L),
      start = gs, end = ge, pam_start = pstart, stringsAsFactors = FALSE)
  }
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(strand = character(), protospacer = character(), pam = character(),
               start = integer(), end = integer(), pam_start = integer(),
               stringsAsFactors = FALSE)
  res <- res[order(res$start, res$strand), , drop = FALSE]
  rownames(res) <- NULL
  res
}

#' Position of a genomic base within a protospacer
#' @noRd
.pos_in_guide <- function(target, start, end, strand) {
  ifelse(strand == "+", target - start + 1L, end - target + 1L)
}
