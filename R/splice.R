#' Locate splice sites in a record
#'
#' Returns the two invariant intronic bases of every usable splice site: the
#' `GT` opening each intron (donor) and the `AG` closing it (acceptor). The
#' acceptor of the first exon and the donor of the last exon are not emitted,
#' since neither corresponds to a real junction.
#'
#' @param record An `hla_record`.
#' @return A data frame with columns `exon`, `site` (`"donor"` or
#'   `"acceptor"`), `base1`, `base2` (genomic positions of the two invariant
#'   bases) and `dinucleotide`.
#' @examples
#' splice_sites(example_record())
#' @export
splice_sites <- function(record) {
  ex <- exons(record)
  n  <- nrow(ex)
  if (n < 2L) return(data.frame(exon = integer(), site = character(),
                                base1 = integer(), base2 = integer(),
                                dinucleotide = character(), stringsAsFactors = FALSE))
  don <- data.frame(exon = ex$index[-n], site = "donor",
                    base1 = ex$end[-n] + 1L, base2 = ex$end[-n] + 2L,
                    stringsAsFactors = FALSE)
  acc <- data.frame(exon = ex$index[-1L], site = "acceptor",
                    base1 = ex$start[-1L] - 2L, base2 = ex$start[-1L] - 1L,
                    stringsAsFactors = FALSE)
  out <- rbind(don, acc)
  out <- out[out$base1 >= 1L & out$base2 <= nchar(record$genomic), , drop = FALSE]
  out$dinucleotide <- substring(record$genomic, out$base1, out$base2)
  out <- out[order(out$base1), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Guides that place a splice-site base in the editing window
#'
#' Enumerates every protospacer in the gene and keeps those covering one of the
#' invariant intronic bases of a splice site. One row is returned per
#' protospacer and per covered base. `target_base` is that base as read on the
#' protospacer's own strand, which determines the chemistry: a `C` is editable
#' by a cytosine base editor, an `A` by an adenine base editor.
#'
#' @param record An `hla_record`.
#' @param pam IUPAC PAM motif.
#' @param guide_length Integer protospacer length.
#' @return A data frame with columns `exon`, `site`, `target`, `strand`,
#'   `protospacer`, `pam`, `start`, `end`, `position`, `target_base` and
#'   `editor`.
#' @examples
#' head(splice_guides(example_record()))
#' @export
splice_guides <- function(record, pam = "NGG", guide_length = 20L) {
  stopifnot(inherits(record, "hla_record"))
  ss <- splice_sites(record)
  g  <- find_protospacers(record$genomic, pam = pam, guide_length = guide_length)
  empty <- data.frame(exon = integer(), site = character(), target = integer(),
                      strand = character(), protospacer = character(),
                      pam = character(), start = integer(), end = integer(),
                      position = integer(), target_base = character(),
                      editor = character(), stringsAsFactors = FALSE)
  if (!nrow(ss) || !nrow(g)) return(empty)

  tg <- rbind(
    data.frame(exon = ss$exon, site = ss$site, target = ss$base1, stringsAsFactors = FALSE),
    data.frame(exon = ss$exon, site = ss$site, target = ss$base2, stringsAsFactors = FALSE))
  idx <- expand.grid(gi = seq_len(nrow(g)), ti = seq_len(nrow(tg)))
  hit <- tg$target[idx$ti] >= g$start[idx$gi] & tg$target[idx$ti] <= g$end[idx$gi]
  idx <- idx[hit, , drop = FALSE]
  if (!nrow(idx)) return(empty)

  out <- data.frame(
    exon = tg$exon[idx$ti], site = tg$site[idx$ti], target = tg$target[idx$ti],
    strand = g$strand[idx$gi], protospacer = g$protospacer[idx$gi],
    pam = g$pam[idx$gi], start = g$start[idx$gi], end = g$end[idx$gi],
    stringsAsFactors = FALSE)
  out$position <- .pos_in_guide(out$target, out$start, out$end, out$strand)
  out$target_base <- substring(out$protospacer, out$position, out$position)
  out$editor <- c(C = "CBE", A = "ABE")[out$target_base]
  out$editor[is.na(out$editor)] <- NA_character_
  out <- out[order(out$target, out$strand, out$position), , drop = FALSE]
  rownames(out) <- NULL
  out
}
