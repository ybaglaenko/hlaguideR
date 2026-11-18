test_that("cds_to_genomic is base-exact across the whole CDS", {
  rec <- example_record()
  p <- seq_len(nchar(rec$coding))
  gpos <- cds_to_genomic(p, rec)
  expect_false(anyNA(gpos))
  expect_equal(substring(rec$coding, p, p),
               substring(rec$genomic, gpos, gpos))
})

test_that("stop-codon guides target a cytosine on the guide strand", {
  g <- stop_codon_guides(example_record())
  expect_true(nrow(g) > 0)
  expect_equal(unique(substring(g$protospacer, g$position, g$position)), "C")
  expect_true(all(g$codon %in% c("CAG", "CGA", "CAA", "TGG")))
  expect_true(all(g$target >= g$start & g$target <= g$end))
})

test_that("NMD escape follows the last-junction rule", {
  rec <- example_record()
  g <- stop_codon_guides(rec, nmd_rule_nt = 55L)
  expect_equal(g$nmd_escape, g$distance_to_last_junction < 55L)
})
