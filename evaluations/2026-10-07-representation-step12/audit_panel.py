#!/usr/bin/env python3
"""Classify current panel deficits without changing production assignments."""
import argparse
import csv
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--contract', type=Path, default=ROOT/'test_data/tmp_representation_step12/accepted-final/windows/0/gap-contract.tsv')
parser.add_argument('--pgphase-output', type=Path, default=ROOT/'test_data/tmp_representation_step12/accepted-final/full-frozen')
args = parser.parse_args()
spec = importlib.util.spec_from_file_location('read_audit', ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
pg_tags, pg_status = audit.assignments(args.pgphase_output/'phased.bam', truth)
hi_tags, hi_status = json.loads((ROOT/'test_data/tmp_gap_fix78/hi-tags.json').read_text())
with args.contract.open() as handle:
    contract = {}
    for row in csv.DictReader(handle, delimiter='\t'):
        if row['window'] in contract:
            assert contract[row['window']] == row
        contract[row['window']] = row
with (ROOT/'evaluations/2026-09-16-test-panel/panel.tsv').open() as handle:
    panel = list(csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t'))
assert {f"{r['gap_left']}-{r['gap_right']}" for r in panel} == set(contract)

def classify(spans, scorable, correct, core, hi_correct, hi_core):
    if not spans:
        return 'no_spanning_block'
    if correct * 5 < scorable * 4:
        return 'below_80pct'
    if correct < hi_correct:
        return 'total_correct_deficit'
    if core < hi_core:
        return 'core_only_deficit'
    return 'passes'

def extents(path):
    result = {}
    with pysam.VariantFile(str(path)) as vcf:
        for row in vcf:
            sample = next(iter(row.samples.values()))
            gt, ps = sample.get('GT'), sample.get('PS')
            if not sample.phased or not gt or None in gt or len(set(gt)) < 2 or not ps:
                continue
            alt = row.alts[0]
            shared = 0
            while shared < min(len(row.ref), len(alt)) and row.ref[shared] == alt[shared]:
                shared += 1
            left, right = result.get(ps, (row.pos, row.pos + shared))
            result[ps] = min(left, row.pos), max(right, row.pos + shared)
    return result

blocks = extents(args.pgphase_output/'phased.vcf')
candidates = list(csv.DictReader((args.pgphase_output/'candidates.tsv').open(), delimiter='\t'))
results = []
members = set()
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as original:
    for row in panel:
        left, right = int(row['gap_left']), int(row['gap_right'])
        key = f'{left}-{right}'
        names = {r.query_name for r in original.fetch('CHM13#0#chr20', left-1, right)
                 if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
        members.update(names)
        tools = {}
        for tool, tags, status in [('pgphase', pg_tags, pg_status), ('hiphase', hi_tags, hi_status)]:
            counts = Counter(status.get(q, 'unphased') for q in names)
            cores = Counter(tags[q][1] for q in names if status.get(q) == 'correct' and 0 < tags[q][1] < 1_000_000_000)
            tools[tool] = {'counts': dict(counts), 'core_correct': max(cores.values(), default=0)}
        native = contract[key]
        assert int(native['scorable']) == len(names)
        assert int(native['hiphase_correct']) == tools['hiphase']['counts'].get('correct', 0)
        assert int(native['hiphase_core_correct']) == tools['hiphase']['core_correct']
        category = classify(int(native['spans']), len(names), int(native['correct']), int(native['core_correct']),
                            int(native['hiphase_correct']), int(native['hiphase_core_correct']))
        spans = any(lo <= left and hi >= right for lo, hi in blocks.values())
        pg, hi = tools['pgphase'], tools['hiphase']
        full_category = classify(spans, len(names), pg['counts'].get('correct', 0), pg['core_correct'],
                                 hi['counts'].get('correct', 0), hi['core_correct'])
        eligible = hi['counts'].get('correct', 0) * 5 >= len(names) * 4
        sites = [r for r in candidates if left < int(r['POS']) < right]
        rescued = sum(pg_status.get(q) == 'correct' and pg_tags[q][1] >= 1_000_000_000 for q in names)
        missed = Counter(pg_status.get(q, 'unphased') for q in names if hi_status.get(q) == 'correct')
        results.append({'window': key, 'scorable': len(names), 'hiphase_80pct': eligible,
                        'native_category': category, 'native': native, 'full_category': full_category,
                        'full_spans': spans, 'tools': tools, 'correct_rescues': rescued,
                        'pgphase_status_of_hiphase_correct_reads': dict(missed),
                        'sites_inside_gap': len(sites), 'site_categories': dict(Counter(r['CATEGORY'] for r in sites)),
                        'unphased_sites_inside_gap': sum(int(r['PHASE_SET']) <= 0 for r in sites)})

summary = {'panel_windows': len(results),
           'native_categories': dict(Counter(r['native_category'] for r in results)),
           'native_hiphase_80pct_categories': dict(Counter(r['native_category'] for r in results if r['hiphase_80pct'])),
           'full_categories': dict(Counter(r['full_category'] for r in results)),
           'full_hiphase_80pct_categories': dict(Counter(r['full_category'] for r in results if r['hiphase_80pct'])),
           'native_vs_full_category_changes': sum(r['native_category'] != r['full_category'] for r in results),
           'unique_original_primary_truth_reads_in_panel': len(members)}
eligible_open = [r for r in results if r['hiphase_80pct'] and r['full_category'] != 'passes']
summary['eligible_open_interior_evidence'] = {
    'no_sites': sum(not r['sites_inside_gap'] for r in eligible_open),
    'all_sites_unphased': sum(r['sites_inside_gap'] > 0 and r['unphased_sites_inside_gap'] == r['sites_inside_gap'] for r in eligible_open),
    'repeat_or_noisy_sites_only': sum(r['sites_inside_gap'] > 0 and set(r['site_categories']) <= {'REP_HET_INDEL', 'NOISY_CAND_HET'} for r in eligible_open)}
report = {'production_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(),
          'pgphase_output': str(args.pgphase_output), 'contract': str(args.contract),
          'hiphase_tags': 'test_data/tmp_gap_fix78/hi-tags.json (existing globally oriented cache)',
          'summary': summary, 'windows': results}
(OUT/'panel-audit.json').write_text(json.dumps(report, indent=2)+'\n')
(OUT/'native-gap-contract.tsv').write_text(args.contract.read_text())
print(json.dumps(summary, indent=2))
