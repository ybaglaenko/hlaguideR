#' Cas12a protospacers
#'
#' Enumerates protospacers downstream of a 5' PAM across both strands. No
#' editing window applies, since Cas12a is used here for cleavage or detection
#' rather than base editing.
#'
#' @param record An `hla_record`.
#' @param pam IUPAC PAM motif; `"TTTV"` by default.
#' @param guide_length Integer protospacer length.
#' @return A data frame as returned by [find_protospacers()].
#' @examples
#' cas12a_guides(example_record())
#' @export
cas12a_guides <- function(record, pam = "TTTV", guide_length = 20L) {
  stopifnot(inherits(record, "hla_record"))
  find_protospacers(record$genomic, pam = pam, guide_length = guide_length,
                    pam_side = "5prime")
}

#' Introduce deliberate seed mismatches
#'
#' Substitutes each seed-region base with its complement, one position at a
#' time, producing one variant per position. Mismatches in the seed are poorly
#' tolerated, so these variants are used to sharpen discrimination between
#' closely related alleles.
#'
#' @param protospacer Character vector of protospacers.
#' @param positions Integer vector of positions to mutate.
#' @return A data frame with columns `protospacer` (the original), `position`
#'   and `variant`.
#' @examples
#' seed_variants("ACGTACGTACGTACGTACGT", positions = 1:3)
#' @export
seed_variants <- function(protospacer, positions = 1:5) {
  if (!is.character(protospacer))
    stop("`protospacer` must be a character vector.", call. = FALSE)
  positions <- as.integer(positions)
  do.call(rbind, lapply(protospacer, function(p) {
    pos <- positions[positions >= 1L & positions <= nchar(p)]
    if (!length(pos)) return(NULL)
    v <- vapply(pos, function(i) {
      ch <- strsplit(p, NULL)[[1]]
      ch[i] <- chartr("ACGT", "TGCA", ch[i])
      paste(ch, collapse = "")
    }, character(1))
    data.frame(protospacer = p, position = pos, variant = v,
               stringsAsFactors = FALSE)
  }))
}
