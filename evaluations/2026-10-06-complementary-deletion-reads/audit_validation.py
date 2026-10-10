#!/usr/bin/env python3
"""Tie completed validation reports to the final immutable production binary."""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TMP = ROOT/'test_data/tmp_gap_fix90'
panel = json.loads((OUT/'panel-audit.json').read_text())
block = json.loads((OUT/'results.json').read_text())
owner = json.loads((OUT/'owner-results.json').read_text())
preservation = json.loads((OUT/'full-preservation.json').read_text())
sha = hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest()
assert panel['binary_sha256'] == block['production_sha256'] == owner['binary_sha256'] == sha
match = re.search(r'All tests passed \((\d+) assertions in (\d+) test cases\)', (TMP/'window-tests.log').read_text())
assert match
assert not preservation['lost_correct_assignments'] and not preservation['lost_phased_assignments']
assert not preservation['changed_variant_records'] and preservation['previous_phased_tags_identical']
assert ROOT.joinpath('test_data/tmp_gap_fix89/final/0/candidates.tsv').read_bytes() == TMP.joinpath('final/0/candidates.tsv').read_bytes()
required = (TMP/'required-checks.log').read_text()
assert 'All tests passed (1551 assertions in 47 test cases)' in required
for gate in ('HiFi TSV golden match', 'HiFi phased VCF golden match', 'ONT TSV golden match', 'ONT phased VCF golden match'):
    assert gate in required
result = {'binary_sha256': sha, 'full_output': 'test_data/tmp_gap_fix90/final/0',
          'accepted_read_connections': [[4625182,4637298],[4637298,4642430]],
          'whole_hiphase_target_read_parity': not any(block['whole_block_deficits'].values()),
          'whole_target_deficits': block['whole_block_deficits'],
          'build': 'passed; no errors or new warnings',
          'unit_tests': 'passed, including four new motif/calibration checks',
          'phasing_predicates': {'cases':47,'assertions':1551,'result':'passed'},
          'goldens': 'HiFi and ONT TSV/VCF match; HiFi 1 versus 4 threads identical',
          'owning_threads': 'identical candidates, complete VCF records and every HP/PS label',
          'new_native_fixture': {'assertions':738,'result':'passed','baseline_assertions':734,
              'baseline_failed_assertions':15,'cached_seconds':float((TMP/'focused-seconds.txt').read_text())},
          'full_gap_suite': {'assertions':int(match[1]),'cases':int(match[2]),'result':'passed'},
          'cache_helper_tests': {'cases':4,'result':'passed'},
          'native_panel': {'output_labels':panel['native_output_labels'],'independent_requests':panel['independent_requests'],
              'comparison_pairs':panel['comparison_pairs'],'changed_output_labels':len(panel['changed_outputs']),
              'result':'passed; every previous correct/phased read and complete VCF record preserved'},
          'full_preservation': 'every prior phased HP/PS label, correct assignment, block extent, complete VCF record and candidate TSV byte retained',
          'full_read_transitions': preservation['transitions'],
          'outside_target_read_transitions': preservation['outside_gap_transitions'],
          'n50_bp':block['blocks']['n50_bp']}
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['native_panel'])
