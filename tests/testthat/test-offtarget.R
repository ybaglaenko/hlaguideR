test_that("Cas-OFFinder input has the documented layout", {
  f <- tempfile()
  write_casoffinder_input(c("ACGTACGTACGTACGTACGT"), "panel.fa", f,
                          pam = "NGG", mismatches = 3L)
  l <- readLines(f)
  expect_equal(l[1], "panel.fa")
  expect_equal(l[2], paste0(strrep("N", 20), "NGG"))
  expect_equal(l[3], "ACGTACGTACGTACGTACGTNNN 3")
})

test_that("output parsing is case-insensitive for the protospacer", {
  f <- tempfile()
  writeLines(c(
    "ACGTACGTACGTACGTACGTNNN\tA*01:01\t10\tACGTACGTACGTACGTACGTAGG\t+\t0",
    "ACGTACGTACGTACGTACGTNNN\tB*07:02\t20\tACGTACGTACGTACGTAcGTAGG\t+\t1"), f)
  d <- read_casoffinder(f)
  expect_equal(unique(d$protospacer), "ACGTACGTACGTACGTACGT")
  expect_equal(d$mismatch, c(0L, 1L))
  expect_equal(d$mismatch_positions[[1]], integer())
  expect_equal(d$mismatch_positions[[2]], 18L)
})

test_that("summaries count distinct alleles and scale by the universe", {
  hits <- data.frame(
    protospacer = "AAA",
    allele = c("A*01:01", "A*01:01", "A*02:01", "B*07:02"),
    stringsAsFactors = FALSE)
  uni <- c("A*01:01", "A*02:01", "A*03:01", "B*07:02")
  s <- offtarget_summary(hits, uni)
  expect_equal(s$n_alleles[s$locus == "A"], 2L)
  expect_equal(s$pct_universe[s$locus == "A"], 200 / 3)
  tot <- offtarget_totals(hits, uni)
  expect_equal(tot$n_alleles, 3L)
  expect_equal(tot$pct_universe, 75)
})

test_that("preserving_guides removes every row for an offending protospacer", {
  guides <- data.frame(protospacer = c("AAA", "CCC"), stringsAsFactors = FALSE)
  hits <- data.frame(protospacer = "AAA", allele = "B*07:02",
                     stringsAsFactors = FALSE)
  expect_equal(preserving_guides(guides, hits, "^B\\*07:02")$protospacer, "CCC")
})
