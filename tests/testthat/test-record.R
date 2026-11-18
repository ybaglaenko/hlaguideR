test_that("example record is internally consistent", {
  rec <- example_record()
  expect_s3_class(rec, "hla_record")
  expect_equal(nchar(rec$genomic), 420L)
  expect_equal(nchar(rec$coding), 210L)
  ex <- exons(rec)
  expect_equal(nrow(ex), 3L)
  expect_equal(ex$end, ex$start + ex$length - 1L)
  expect_equal(sum(ex$length), nchar(rec$coding))
})

test_that("exon concatenation reproduces the coding sequence", {
  rec <- example_record()
  ex <- exons(rec)
  expect_equal(paste(substring(rec$genomic, ex$start, ex$end), collapse = ""),
               rec$coding)
})

test_that("cDNA-only API records are rejected with a clear message", {
  x <- list(sequence = list(genomic = "", coding = "ACG"),
            allele_history = list(name = "A*68:126"))
  expect_error(as_hla_record(x), "no genomic sequence")
})

test_that("bad input is caught", {
  expect_error(hla_record("A", "", data.frame()), "non-empty")
  expect_error(hla_record("A", "ACGT", data.frame(x = 1)), "columns")
})
