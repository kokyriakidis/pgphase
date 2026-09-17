# Is the graph-site loading robust and accurate?

Audited before changing anything, because a redesign of a loader that is already
correct trades a known state for an unknown one.

## Accuracy: measured, and it is correct on this catalog

Over `chr20:5,259,406-5,395,085`, the region the target window's run loads:

| property | count |
|---|---:|
| records | 2,195 |
| positions carrying more than one record | **0** |
| records with a `*` spanning-deletion allele | 0 |
| records with a symbolic `<...>` allele | 0 |
| records with characters outside ACGTN | 0 |
| records with lower-case REF or ALT | 0 |
| out-of-order POS transitions | 0 |

Every record is parsed and every site is eligible: the run reports
`2195 data lines, 2195 sites parsed, 2195 eligible`. Injection also honours the
verdict -- `hybrid_inject.cpp:349` skips a site whose structural validation
failed -- so an ineligible site cannot leak into the candidate table.

**Efficiency is not the constraint.** The loader is tabix-region-filtered
(`load_sites_for_region` builds a single `RegionFilter` and iterates only that
span), so a chunk reads only the records it needs: 2,195 for this window,
against a run that takes tens of seconds. Streaming the whole file happens only
when no region filter is set.

## Robustness: four holes, none of which this catalog exercises

Being correct on one file is not the same as being robust, and every hole below
is a silent or fatal failure on a file we have not been handed yet.

| hole | what happened before | now |
|---|---|---|
| `std::stoll` on POS and END | a malformed value **throws out of the parser** and aborts the whole run, naming no line | the record is skipped and counted in `bad_position` |
| symbolic alleles (`<DEL>`, `<INS>`) | passed through to `vcf_to_variant_key`, which treats the angle brackets as sequence and mints a candidate whose ALT is literally `<DEL>` | the record is skipped and counted in `unsupported_allele` |
| lower-case alleles | a soft-masked region writes `acgt`; the allele would never match a candidate derived from upper-case sequence | REF and ALT are upper-cased at load |
| silent discards | a short line returned with no record and no trace, so a catalog that lost half its sites looked exactly like one that loaded cleanly | every discard is counted, and `GraphSiteLoadStats::summary()` is printed per chunk at `--verbose 1` |

`GraphSiteLoadStats` carries `data_lines`, `sites_parsed`, `short_line`,
`bad_position`, `unsupported_allele` and `by_skip_reason`, so the load is
described by the loader rather than inferred from its downstream effect.

## Verification

`test_graph_sites.cpp` feeds a VCF containing one well-formed record, a POS that
is not a number, a POS of zero, a symbolic allele, a lower-case record, a short
line, and a non-numeric END. It asserts all seven lines are counted, exactly the
two well-formed records survive, the three bad positions and the one symbolic
allele and the one short line are each counted, and the lower-case record is
stored upper-cased. Before this change that file **aborted the run** on its
second line.

Validated by fault injection: removing the allele screen makes the test fail
with `the symbolic allele is rejected`, so the assertion is load-bearing.

Behaviour on real data is unchanged: for the target window the candidate table
and the phased VCF are **byte-identical** to the pre-change run. Unit tests 4/4,
injection tests 127, window tests 178.

## What was deliberately not changed

- **No dedup of records at one position.** This catalog has none, so a dedup
  rule would be written against no evidence about which record to prefer.
- **No string_view parse.** It would cut allocations, but the loader reads only
  the chunk's own records and is nowhere near the cost of a run, so the change
  would buy nothing measurable and risk the parse.
