# Rejected broad normalized observation admission

The provisional binary was
`e6707eeb294269e7f507ac5ebb1384a3e37412230648cd07078357a7664c4c9c`.
It allowed uniquely normalized, covered BAM calls to fill missing graph calls.
All 141 native panel measurements were unchanged. The full chromosome gained
33 output rescues: 29 correct, four discordant, replacing previously unphased
tags. The same 33 reads have HiPhase counts 29 correct, one discordant and three
unphased. All established tags and connected cores remained unchanged.

Nevertheless only 105 of 107 named checks passed. The owning 47 Mb deletion
regression gained one discordant read (37 vs ceiling 36; concordance 0.990736
vs floor 0.9909), and the owning 33 Mb short-insertion regression gained one
(133 vs ceiling 132). Those stronger existing bounds are retained. The added
normalized calls were usable without an independent graph observation; sequence
identity and alignment bounds alone did not certify source reference-class
confidence. No new gap closure was accepted.

The accepted integration restricts normalized observations to reads already
observed independently by the graph and keeps conflicting channel calls visible.
Missing normalized graph calls await complete allele contrasts, read confidence
and molecule/provenance handling. Fast fixtures explicitly gate this distinction.
The provisional full output is retained in
`test_data/tmp_representation_step2/full-provisional/`, the binary as
`pgphase.provisional` in that work directory, and the two failed tests are saved
here. Final verification uses a rebuilt binary and separate replay state; no
expectations are refreshed.
