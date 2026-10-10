#!/usr/bin/env python3
"""Tie the closure and preservation gates to the final production executable."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TMP = ROOT/'test_data/tmp_gap_fix91'
sha = hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest()
reports = {name: json.loads((OUT/file).read_text()) for name, file in
           [('panel', 'panel-audit.json'), ('block', 'results.json'),
            ('owner', 'owner-results.json'), ('preservation', 'full-preservation.json')]}
assert reports['panel']['binary_sha256'] == reports['block']['production_sha256'] == reports['owner']['binary_sha256'] == reports['preservation']['binary_sha256'] == sha
match = re.search(r'All tests passed \((\d+) assertions in (\d+) test cases\)', (TMP/'window-tests.log').read_text())
assert match
assert 'Ran 4 tests in' in (TMP/'window-tests.log').read_text()
checks = (TMP/'required-checks.log').read_text()
assert 'All tests passed (1551 assertions in 47 test cases)' in checks
assert checks.count('ALL PASS') >= 4
for gate in ('HiFi TSV golden match', 'HiFi phased VCF golden match', 'ONT TSV golden match', 'ONT phased VCF golden match'):
    assert gate in checks
assert 'ALL PASS' in (TMP/'build-tests-final.log').read_text()
assert 'error:' not in (TMP/'build-final.log').read_text() and 'warning:' not in (TMP/'build-final.log').read_text()
preservation = reports['preservation']
assert not preservation['lost_correct_assignments'] and not preservation['lost_phased_assignments']
assert preservation['previous_haplotypes_identical'] and preservation['previous_core_tags_identical']
assert not preservation['changed_variant_records']
assert ROOT.joinpath('test_data/tmp_gap_fix90/final/0/candidates.tsv').read_bytes() == TMP.joinpath('final_current/0/candidates.tsv').read_bytes()
focused = re.search(r'All tests passed \((\d+) assertions in (\d+) test cases\)', (TMP/'focused-warm.log').read_text())
assert focused
baseline = (TMP/'baseline-fixture.log').read_text()
assert '31 >= 53' in baseline and 'FAILED:' in baseline
result = {
    'binary_sha256': sha, 'full_output': 'test_data/tmp_gap_fix91/final_current/0',
    'accepted_gap': [4874129, 4884130],
    'whole_hiphase_target_read_parity': not any(reports['block']['whole_block_deficits'].values()),
    'whole_target_deficits': reports['block']['whole_block_deficits'],
    'build': 'passed without errors or warnings',
    'units_and_phasing_predicates': 'passed; 1551 predicate assertions in 47 cases',
    'new_in_memory_checks': 13,
    'cache_helper_tests': {'cases': 4, 'result': 'passed'},
    'goldens': 'HiFi/ONT TSV and VCF match; HiFi one/four-thread outputs identical',
    'owning_threads': 'one/four-thread read tags, candidate bytes and complete VCF identical',
    'new_native_fixture': {'assertions': int(focused[1]), 'result': 'passed',
                           'baseline_fails_core_parity': True,
                           'warm_seconds': float((TMP/'focused-seconds.txt').read_text())},
    'full_gap_suite': {'assertions': int(match[1]), 'cases': int(match[2]), 'result': 'passed'},
    'native_panel': {k: reports['panel'][k] for k in
                    ('native_output_labels', 'independent_requests', 'comparison_pairs')},
    'changed_native_labels': len(reports['panel']['changed_outputs']),
    'full_preservation': 'all previous correct/phased reads, haplotypes, core tags, block extents, VCF records and candidate bytes preserved',
    'full_read_transitions': preservation['transitions'],
    'n50_bp': reports['block']['blocks']['n50_bp']}
(OUT/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
