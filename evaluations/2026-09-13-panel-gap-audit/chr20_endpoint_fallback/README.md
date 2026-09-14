# chr20 endpoint-scoped phase fallback

This replay uses the 114-case `chr20_clean_scoped_iterations` manifest and
changes only the projection of each case's two audited endpoint positions.
The projector first uses ordinary HP/PS support, then an exact phased native
pgphase call, then a phased input-caller block with at least three heterozygous
records. Caller blocks use a separate PS namespace and cannot merge pgphase
blocks.

The replay changes two records. The native addition at 26,602,087 agrees with
the established orientation of its pgphase block against assembly truth. The
caller addition at 26,626,247 is isolated from pgphase blocks. Target status is
88 joined, 26 split, and zero endpoint-unphased. The shared 19 kb phase gap is
still open.

The production panel command is recorded in `run_link_panel.py`; it supplies
`--native-vcf`, `--retain-caller-phase-blocks`, and the two `--fallback-site`
positions from each manifest row.
