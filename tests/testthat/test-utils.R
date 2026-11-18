test_that("revcomp round-trips and handles IUPAC", {
  expect_equal(revcomp("ACGT"), "ACGT")
  expect_equal(revcomp(revcomp("ACGTTGCAN")), "ACGTTGCAN")
  expect_equal(revcomp("TTTV"), "BAAA")
})

test_that("iupac_regex expands ambiguity codes", {
  expect_equal(iupac_regex("NGG"), "[ACGT]GG")
  expect_equal(iupac_regex("TTTV"), "TTT[ACG]")
  expect_error(iupac_regex("XYZ"), "IUPAC")
})

test_that("protospacer coordinates are plus-strand and consistent", {
  rec <- example_record()
  g <- find_protospacers(rec$genomic, "NGG", 20L)
  expect_true(nrow(g) > 0)
  expect_true(all(g$end - g$start + 1L == 20L))
  expect_true(all(g$start >= 1L & g$end <= nchar(rec$genomic)))
  span <- substring(rec$genomic, g$start, g$end)
  expect_equal(ifelse(g$strand == "+", span, revcomp(span)), g$protospacer)
})

test_that("overlapping PAMs are all reported", {
  # a terminal GGG yields two NGG matches (AGG, GGG); GGGG yields three
  g2 <- find_protospacers(paste0(strrep("A", 22), "GGG"), "NGG", 20L, strands = "+")
  expect_equal(nrow(g2), 2L)
  g3 <- find_protospacers(paste0(strrep("A", 21), "GGGG"), "NGG", 20L, strands = "+")
  expect_equal(nrow(g3), 3L)
  expect_equal(g3$start, 1:3)
})

test_that("guides running off the sequence end are dropped", {
  # PAM present but fewer than guide_length bases upstream
  g <- find_protospacers(paste0(strrep("A", 5), "AGG"), "NGG", 20L, strands = "+")
  expect_equal(nrow(g), 0L)
})

test_that("5-prime PAMs sit upstream of the protospacer", {
  s <- paste0("TTTA", strrep("C", 20))
  g <- find_protospacers(s, "TTTV", 20L, pam_side = "5prime", strands = "+")
  expect_equal(g$protospacer, strrep("C", 20))
  expect_equal(g$pam, "TTTA")
})
