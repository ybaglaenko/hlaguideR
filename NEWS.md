# hlaguideR 0.1.0

First release.

## Guide design

* `splice_guides()` enumerates protospacers that place an invariant splice-site
  base inside the guide, reporting one row per protospacer and per covered base.
  The base editor is derived from the sequence rather than assumed: `target_base`
  is the splice base as read on the protospacer's own strand, so a `C` is
  editable by a cytosine base editor and an `A` by an adenine base editor.
* `stop_codon_guides()` finds codons convertible to a stop by a single
  C-to-T transition (`CAG`, `CGA`, `CAA`, and `TGG` through its antisense `CCA`)
  and returns the guides that place the required cytosine inside the
  protospacer. Nonsense-mediated decay is reported through
  `distance_to_last_junction` and `nmd_escape`, not silently filtered.
* `cas12a_guides()` enumerates protospacers downstream of a 5' PAM, and
  `seed_variants()` introduces deliberate single-base seed mismatches for
  allele discrimination.
* `find_protospacers()` underlies all three. It scans both strands, reports
  every overlapping PAM occurrence, and always returns plus-strand coordinates
  with `start <= end`; the `strand` column records orientation.

## Off-target accounting

* `write_casoffinder_input()` and `read_casoffinder()` exchange files with
  Cas-OFFinder without depending on the binary, so the package installs and
  checks anywhere. `read_casoffinder()` resolves mismatch positions at parse
  time and returns both the verbatim and an upper-cased sequence, so no
  downstream join depends on letter case.
* `offtarget_summary()` and `offtarget_totals()` count distinct alleles hit per
  protospacer, by locus and in total, optionally scaled by the panel searched.
* `preserving_guides()` drops every protospacer that hits a protected allele.

## Records

* `hla_record()` bundles a genomic sequence with its exon/intron feature table;
  `cds_to_genomic()` maps coding positions onto genomic coordinates.
* `as_hla_record()` converts an IPD-IMGT/HLA API response and rejects
  cDNA-only records with an explicit message, rather than failing later in the
  pipeline. Many alleles in that database carry no genomic sequence.
* `example_record()` returns a small synthetic gene used by the examples and
  tests, so neither requires network access or redistributed database sequence.
