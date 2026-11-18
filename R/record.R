#' Create an HLA gene record
#'
#' An `hla_record` bundles a genomic sequence with the exon/intron feature
#' table that indexes it. All guide-design functions take a record.
#'
#' @param allele Character scalar naming the allele, e.g. `"B*27:05"`.
#' @param genomic Character scalar; the full genomic sequence.
#' @param features A data frame with columns `type` (one of `"exon"`,
#'   `"intron"`, `"5utr"`, `"3utr"`), `number`, `start` (1-based position in
#'   `genomic`) and `length`.
#' @param coding Optional coding sequence. When `NULL` it is built by
#'   concatenating the exons.
#' @param codon_start Integer 1, 2 or 3; the first complete codon's offset
#'   within `coding`.
#' @return An object of class `hla_record`.
#' @examples
#' rec <- example_record()
#' rec
#' @export
hla_record <- function(allele, genomic, features, coding = NULL, codon_start = 1L) {
  if (!is.character(allele) || length(allele) != 1L)
    stop("`allele` must be a character scalar.", call. = FALSE)
  if (!is.character(genomic) || length(genomic) != 1L || !nzchar(genomic))
    stop("`genomic` must be a non-empty character scalar.", call. = FALSE)
  need <- c("type", "number", "start", "length")
  if (!is.data.frame(features) || !all(need %in% names(features)))
    stop("`features` needs columns: ", paste(need, collapse = ", "), call. = FALSE)
  codon_start <- as.integer(codon_start)
  if (!codon_start %in% 1:3) stop("`codon_start` must be 1, 2 or 3.", call. = FALSE)

  genomic <- toupper(genomic)
  features$type <- tolower(as.character(features$type))
  features <- features[order(features$start), , drop = FALSE]
  rownames(features) <- NULL
  ex <- features[features$type == "exon", , drop = FALSE]
  if (!nrow(ex)) stop("`features` contains no exon rows.", call. = FALSE)
  if (max(ex$start + ex$length - 1L) > nchar(genomic))
    stop("Exon coordinates run past the end of `genomic`.", call. = FALSE)

  built <- paste(substring(genomic, ex$start, ex$start + ex$length - 1L), collapse = "")
  if (is.null(coding)) {
    coding <- built
  } else {
    coding <- toupper(coding)
    if (!identical(coding, built))
      warning("`coding` differs from the concatenated exons; using `coding` as given.",
              call. = FALSE)
  }
  structure(list(allele = allele, genomic = genomic, coding = coding,
                 codon_start = codon_start, features = features),
            class = "hla_record")
}

#' @param x An `hla_record`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @rdname hla_record
#' @export
print.hla_record <- function(x, ...) {
  ex <- exons(x)
  cat("<hla_record> ", x$allele, "\n", sep = "")
  cat("  genomic: ", nchar(x$genomic), " nt | coding: ", nchar(x$coding),
      " nt | codon_start: ", x$codon_start, "\n", sep = "")
  cat("  exons:   ", nrow(ex), "\n", sep = "")
  invisible(x)
}

#' Exon rows of a record
#'
#' @param record An `hla_record`.
#' @return A data frame of exon features with an added `end` column holding the
#'   last base of each exon (that is, `start + length - 1`).
#' @examples
#' exons(example_record())
#' @export
exons <- function(record) {
  stopifnot(inherits(record, "hla_record"))
  ex <- record$features[record$features$type == "exon", , drop = FALSE]
  ex <- ex[order(ex$start), , drop = FALSE]
  ex$end <- ex$start + ex$length - 1L
  ex$index <- seq_len(nrow(ex))
  rownames(ex) <- NULL
  ex
}

#' A small synthetic record for examples and tests
#'
#' A 420 nt three-exon gene with canonical `GT`/`AG` intron termini, a 210 nt
#' coding sequence free of internal stop codons, and four codons targetable by
#' cytosine base editing. The sequence is synthetic and carries no database
#' licence.
#'
#' @return An `hla_record`.
#' @examples
#' example_record()
#' @export
example_record <- function() {
  gene <- paste0(
    "TTGCGTGCGCAACGGCCTGTTCAATCGTGGATGGCCACGTTCCTGGAAGACCTCATTCAGGTTACATGCGTGCGCCCTGGACTCTTACGTGTCGCACGGCACACTCACGTCCATTTTACCTGGGATTGCGATAATGGCCACTGGATGCGTTATACTCAGTATTACTCTAGTCCGTGAGTACGCGATCATACCGCAGCCCAGTAATGCTCG",
    "GCGAACGGATGGTTCATCAAACCCACCGGTTTACATATTGTTGCTTGTGTGTTCAACACCTAATCTGAAGAAGTGCATATAATTCGTCGGTCTCGTAGTTACGACGCCAACATGCATGAGCGCGGATATCCTTGGACGTCGCCGTTAAGCGGATGCCTAACTGTGGCTAGTGGCGCATGACGTTGCACCACCCGACCTCTATGCCAGGAC")
  feat <- data.frame(
    type   = c("5utr", "exon", "intron", "exon", "intron", "exon", "3utr"),
    number = c(NA, 1L, 1L, 2L, 2L, 3L, NA),
    start  = c(1L, 31L, 91L, 171L, 261L, 331L, 391L),
    length = c(30L, 60L, 80L, 90L, 70L, 60L, 30L),
    stringsAsFactors = FALSE)
  hla_record("EXAMPLE*01:01", gene, feat)
}

#' Build a record from an IPD-IMGT/HLA API response
#'
#' Converts the list returned by the IPD-IMGT/HLA allele endpoint into an
#' [hla_record()]. Records that carry no genomic sequence -- a large share of
#' the database, since many alleles are typed from cDNA only -- are rejected
#' with an informative error rather than failing later.
#'
#' @param x A list parsed from the API's JSON response.
#' @return An `hla_record`.
#' @export
as_hla_record <- function(x) {
  if (!is.list(x) || is.null(x$sequence))
    stop("Response carries no `sequence` element.", call. = FALSE)
  gen <- x$sequence$genomic
  nm  <- if (!is.null(x$allele_history$name))
    utils::tail(x$allele_history$name, 1L) else NA_character_
  if (is.null(gen) || !nzchar(gen))
    stop("Allele ", nm, " has no genomic sequence (cDNA-only record); ",
         "splice-site design is not possible for it.", call. = FALSE)
  feat <- x$feature$genomic
  if (is.null(feat)) stop("Response carries no genomic feature table.", call. = FALSE)
  feat <- data.frame(type = feat$type, number = feat$number,
                     start = as.integer(feat$start),
                     length = as.integer(feat$length), stringsAsFactors = FALSE)
  cs <- if (is.null(x$sequence$codon_start)) 1L else as.integer(x$sequence$codon_start)
  hla_record(if (is.na(nm)) "unknown" else nm, gen, feat,
             coding = x$sequence$coding, codon_start = cs)
}

#' Fetch an allele record from IPD-IMGT/HLA
#'
#' Requires network access and the \pkg{jsonlite} package.
#'
#' @param accession IPD accession, e.g. `"HLA00001"`.
#' @return An `hla_record`.
#' @examples
#' \donttest{
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   try(fetch_hla_record("HLA00001"))
#' }
#' }
#' @export
fetch_hla_record <- function(accession) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for fetch_hla_record().", call. = FALSE)
  if (!is.character(accession) || length(accession) != 1L)
    stop("`accession` must be a character scalar.", call. = FALSE)
  url <- paste0("https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/", accession)
  res <- tryCatch(jsonlite::fromJSON(url),
                  error = function(e)
                    stop("Could not retrieve ", accession, ": ",
                         conditionMessage(e), call. = FALSE))
  as_hla_record(res)
}
