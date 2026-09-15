# Recovery error attribution

Compared current binary with recovery disabled (`/tmp/pgphase-error-audit-no-recovery/command.json`) against the accepted recovery run (`/tmp/pgphase-confirm-original/command.json`). Same graph/BAM inputs and parental read truth; no DeepVariant input to phasing.

| Run | Evaluated reads | Discordant reads |
| --- | ---: | ---: |
| Recovery disabled | 185835 | 586 |
| Recovery enabled | 189042 | 2427 |

The cohorts differ. `block_changes.tsv` separately compares common reads grouped by original and final PS. Every compared original block preserves a uniform HP-label transformation. Nevertheless, original blocks become reversed relative to the majority of their merged component: PS52434849 goes from 1/489 to 488/489 discordant; PS46663811 goes from 1/235 to 234/235; PS46748667 goes from 0/347 to 347/347. These are real relative-orientation failures, not merely newly assigned noisy reads.

The frozen audit at gap 52394827–52434849 chooses no flip at tier 3, with left [44,0,0,50], right [53,0,0,66]. The original blocks have opposite truth orientations, so a flip is required. The BAM-only view also incorrectly chooses no flip. This demonstrates that excellent flank attachment and agreement between correlated proposal views cannot certify the internal gap bridge. Gap 46844468–46895038 likewise accepts the incorrect no-flip parity, attaching two upstream blocks to a larger component in the wrong orientation. This investigation has not yet fixed those production decisions.

The evaluator additionally mixed maternal and paternal assembly positions when ordering reads. Correct common input-reference ordering changes the recovery run's read transition diagnostics from 351 switches / 414 flips to 205 / 267, without changing discordance. These are read diagnostics, not variant switch/flip metrics. The evaluator now records metric unit, coordinate system and version, sorts ties by read name, separates contigs for transitions and excludes unpositioned reads from transitions. Legacy JSON count keys remain for consumers. Other legacy span metrics still use truth coordinates and require a separate audit before interpreting as reference block contiguity.

Regression `scripts/test_evaluate_phase_accuracy.py` shifts paternal truth positions by 10 kb while holding input read positions and truth labels fixed, asserting unchanged transition diagnostics and discordance.
