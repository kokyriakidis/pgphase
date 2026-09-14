# chr20 46.4 Mb: an available bridge discarded by nearest-link selection

The shared DeepVariant endpoint SNP 46384021 and insertion 46405048 lie in
one correctly oriented competitor block but in different pgphase blocks.
The native clean boundary is 46384021–46407462.

The insertion is present in the graph catalog: its 46405052 AA→AAA allele
left-normalizes to 46405048 T→TA (see catalog.normalized.vcf). The SNP lies
inside a larger catalog snarl and is already a clean BAM candidate. Tier 3
retrieves the MSA-verified insertion at internal key 46405049, with 65 allele
observations (42 reference, 23 alternate). Retrieval is therefore sufficient.

Three independent primary MAPQ60 molecules span the SNP and insertion:

| Read | SNP allele | Insertion allele |
|---|---:|---:|
| m84031_231217_034919_s2/170592309/ccs | 0 | 0 |
| m84031_231217_034919_s2/84544880/ccs | 1 | 1 |
| m84031_231217_034919_s2/40244203/ccs | 0 | 0 |

The same observations occur in the recovery matrix (columns 50 and 65).
Nevertheless, choosing only the nearest sufficient preceding link attaches
this insertion to the right component and discards its link to the left.
A later vertex cannot join two earlier components under that rule.

The recovery-round change retains supported edges within the existing het
window, orders them by net evidence, and uses parity-aware component unions.
It keeps the first-round k-means and the existing block-stitching path.
The reproduced target improves from two blocks to one, with the same 420
truth-evaluated reads and zero discordance or read switchflips. Unit tests
cover the two-component bridge, orientation composition, convergence, and
weaker contradictory cycles. The first local reproduction used a stricter
net-support gate; the final candidate retains the original majority-support
gate and is evaluated separately in chr20_graph_support2.
