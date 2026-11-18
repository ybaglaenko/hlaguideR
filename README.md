# hlaguideR

Design CRISPR base-editing and Cas12a guides against HLA alleles.

`hlaguideR` enumerates protospacers that disrupt splice sites or install
premature termination codons in HLA genes, annotates the position of every
editable base within each protospacer, and summarises off-target searches
across a panel of alleles.

The package has no hard dependencies beyond base R. Alignment search itself is
delegated to an external tool (Cas-OFFinder, Bae, Park and Kim 2014); the
package writes its input and parses its output, so nothing external is needed
to install or run the package.

## Installation

```r
install.packages("hlaguideR_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(hlaguideR)

# a synthetic record, no network required
rec <- example_record()
splice_guides(rec)
stop_codon_guides(rec)

# or a real allele from IPD-IMGT/HLA
rec <- fetch_hla_record("HLA00001")       # A*01:01:01:01
g <- splice_guides(rec)
w <- filter_editing_window(g, min = 3, max = 10)
subset(w, !is.na(editor))                  # guides a base editor can act on
```

`target_base` is the splice-site base as read on the protospacer's own strand,
so it determines the chemistry directly: a `C` is editable by a cytosine base
editor, an `A` by an adenine base editor, and anything else is reported with
`editor = NA` rather than dropped.

## Off-target workflow

```r
write_casoffinder_input(g$protospacer, "panel.fa", "input.txt", pam = "NGG")
# run cas-offinder input.txt C output.txt
hits <- read_casoffinder("output.txt")
offtarget_summary(hits, universe = panel_alleles)
preserving_guides(g, hits, "^B\\*07:02")   # drop guides that hit a protected allele
```

## Scope

This package is the design engine. The interactive application built on top of
it -- population frequency filtering, base-editing efficiency prediction,
primer design and plotting -- lives separately, since it depends on external
binaries and reference databases that cannot ship with a CRAN package.

## License

MIT
