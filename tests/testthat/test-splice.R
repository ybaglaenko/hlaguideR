test_that("splice sites are the canonical GT/AG dinucleotides", {
  ss <- splice_sites(example_record())
  expect_equal(nrow(ss), 4L)
  expect_setequal(ss$dinucleotide[ss$site == "donor"], "GT")
  expect_setequal(ss$dinucleotide[ss$site == "acceptor"], "AG")
})

test_that("first-exon acceptor and last-exon donor are not emitted", {
  ss <- splice_sites(example_record())
  n <- nrow(exons(example_record()))
  expect_false(any(ss$site == "acceptor" & ss$exon == 1L))
  expect_false(any(ss$site == "donor" & ss$exon == n))
})

test_that("splice guides really cover their target base", {
  g <- splice_guides(example_record())
  expect_true(nrow(g) > 0)
  expect_true(all(g$target >= g$start & g$target <= g$end))
  expect_true(all(g$position >= 1L & g$position <= nchar(g$protospacer)))
  expect_equal(substring(g$protospacer, g$position, g$position), g$target_base)
  expect_true(all(g$editor[g$target_base == "C"] == "CBE"))
  expect_true(all(g$editor[g$target_base == "A"] == "ABE"))
  expect_true(all(is.na(g$editor[!g$target_base %in% c("A", "C")])))
})
