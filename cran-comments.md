## Submission summary

New submission. `hlaguideR` enumerates CRISPR protospacers that disrupt splice
sites or install premature termination codons in HLA genes, and summarises
off-target searches across a panel of alleles.

The package has no hard dependencies beyond base R. Alignment search itself is
delegated to an external tool (Cas-OFFinder); the package only writes its input
and parses its output, so nothing external is required to install, check or run
the examples and tests.

## Test environments

* local: aarch64-apple-darwin20.0.0, R 4.5.3 (2026-03-11) --
  `R CMD check --as-cran --no-manual`
* TODO before submitting -- win-builder, R-devel: `devtools::check_win_devel()`
* TODO before submitting -- win-builder, R-release: `devtools::check_win_release()`
* TODO before submitting -- R-hub linux/macos: `rhub::rhub_check()`

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE

      Maintainer: 'Yuriy Baglaenko <yuriy.baglaenko@cchmc.org>'
      New submission

  This is a first submission, so the note is expected.

* checking for future file timestamps ... NOTE

      unable to verify current time

  Raised by the local check machine, which cannot reach the time service used
  for this check. It is unrelated to the package contents.

The local check was run with `--no-manual` because LaTeX is not installed on
that machine; the PDF manual builds on win-builder.

## Downstream dependencies

None. This is a new package, so there are no reverse dependencies to check.
