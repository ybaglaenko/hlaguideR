#' Map a coding-sequence position to a genomic coordinate
#'
#' @param pos Integer vector of 1-based positions in the coding sequence.
#' @param record An `hla_record`.
#' @return Integer vector of genomic coordinates; `NA` where `pos` falls
#'   outside the coding sequence.
#' @examples
#' cds_to_genomic(c(1, 61, 151), example_record())
#' @export
cds_to_genomic <- function(pos, record) {
  ex <- exons(record)
  cds_start <- cumsum(c(1L, utils::head(ex$length, -1L)))
  cds_end   <- cumsum(ex$length)
  vapply(as.integer(pos), function(p) {
    i <- which(p >= cds_start & p <= cds_end)
    if (!length(i)) return(NA_integer_)
    ex$start[i[1L]] + (p - cds_start[i[1L]])
  }, integer(1))
}

#' Guides that install a premature termination codon
#'
#' Finds codons convertible to a stop by a single C-to-T transition --
#' `CAG`, `CGA` and `CAA` on the sense strand, and `TGG` through its antisense
#' `CCA` -- then returns the protospacers that place the required cytosine
#' inside the guide. `position` is the cytosine's index within the protospacer.
#'
#' Nonsense-mediated decay is reported, not filtered: `nmd_escape` is `TRUE`
#' where the new stop codon lies downstream of the final exon-exon junction or
#' within `nmd_rule_nt` bases upstream of it, the range in which transcripts
#' are generally not degraded.
#'
#' @param record An `hla_record`.
#' @param pam IUPAC PAM motif.
#' @param guide_length Integer protospacer length.
#' @param nmd_rule_nt Distance threshold for the last-junction rule.
#' @return A data frame with columns `codon`, `codon_index`, `cds_position`,
#'   `target`, `strand`, `protospacer`, `pam`, `start`, `end`, `position`,
#'   `distance_to_last_junction` and `nmd_escape`.
#' @examples
#' head(stop_codon_guides(example_record()))
#' @export
stop_codon_guides <- function(record, pam = "NGG", guide_length = 20L,
                              nmd_rule_nt = 55L) {
  stopifnot(inherits(record, "hla_record"))
  cds <- record$coding
  off <- record$codon_start - 1L
  n   <- (nchar(cds) - off) %/% 3L
  empty <- data.frame(codon = character(), codon_index = integer(),
                      cds_position = integer(), target = integer(),
                      strand = character(), protospacer = character(),
                      pam = character(), start = integer(), end = integer(),
                      position = integer(), distance_to_last_junction = integer(),
                      nmd_escape = logical(), stringsAsFactors = FALSE)
  if (n < 1L) return(empty)

  s <- off + seq(1L, by = 3L, length.out = n)
  cod <- substring(cds, s, s + 2L)
  # offset of the editable base within each codon
  take <- list(CAG = 0L, CGA = 0L, CAA = 0L, TGG = c(1L, 2L))
  keep <- cod %in% names(take)
  if (!any(keep)) return(empty)
  idx <- which(keep)
  cand <- do.call(rbind, lapply(idx, function(i) {
    o <- take[[cod[i]]]
    data.frame(codon = cod[i], codon_index = i, cds_position = s[i] + o,
               stringsAsFactors = FALSE)
  }))
  cand$target <- cds_to_genomic(cand$cds_position, record)
  cand <- cand[!is.na(cand$target), , drop = FALSE]
  if (!nrow(cand)) return(empty)

  g <- find_protospacers(record$genomic, pam = pam, guide_length = guide_length)
  if (!nrow(g)) return(empty)
  ix <- expand.grid(gi = seq_len(nrow(g)), ci = seq_len(nrow(cand)))
  ok <- cand$target[ix$ci] >= g$start[ix$gi] & cand$target[ix$ci] <= g$end[ix$gi]
  ix <- ix[ok, , drop = FALSE]
  if (!nrow(ix)) return(empty)

  out <- cbind(cand[ix$ci, c("codon", "codon_index", "cds_position", "target")],
               g[ix$gi, c("strand", "protospacer", "pam", "start", "end")])
  out$position <- .pos_in_guide(out$target, out$start, out$end, out$strand)
  base <- substring(out$protospacer, out$position, out$position)
  out <- out[base == "C", , drop = FALSE]
  if (!nrow(out)) return(empty)

  ex <- exons(record)
  last_junction <- if (nrow(ex) > 1L) sum(utils::head(ex$length, -1L)) else NA_integer_
  out$distance_to_last_junction <- if (is.na(last_junction)) NA_integer_ else
    last_junction - out$cds_position
  out$nmd_escape <- !is.na(out$distance_to_last_junction) &
    out$distance_to_last_junction < nmd_rule_nt
  out <- out[order(out$cds_position, out$strand), , drop = FALSE]
  rownames(out) <- NULL
  out
}
