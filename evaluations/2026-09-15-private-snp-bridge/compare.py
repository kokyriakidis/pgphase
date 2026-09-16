#!/usr/bin/env python3
"""Compare two whole-chromosome arms read by read and gap by gap.

Reports the project's regression gate (previously concordant reads that became
discordant) alongside coverage changes and gap-status movement.
"""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True, help='arm directory used as reference')
p.add_argument('--run', type=Path, required=True, help='arm directory under test')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()


def per_read(root):
    with gzip.open(root / 'eval/per_read.tsv.gz', 'rt') as stream:
        return {r['read_name']: r for r in csv.DictReader(stream, delimiter='\t')}


def tiers(root):
    with (root / 'tiers.tsv').open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def gap_status(rows):
    """Best status per gap: a gap is joined if any tier attempt joined it."""
    rank = {'joined': 0, 'split': 1, 'partial': 2, 'vetoed': 3, 'rejected': 4, 'open': 5}
    best = {}
    for row in rows:
        key = (row['CHROM'], row['GAP_LEFT'], row['GAP_RIGHT'])
        status = row['STATUS']
        if key not in best or rank.get(status, 9) < rank.get(best[key], 9):
            best[key] = status
    return best


before, after = per_read(a.baseline), per_read(a.run)
transitions = collections.Counter()
corrupted, recovered, lost = [], [], []
for name, old in before.items():
    new = after.get(name)
    if new is None:
        transitions['read_dropped'] += 1
        lost.append(name)
        continue
    # The evaluator writes 'concordant' lower case and 'DISCORDANT' upper case.
    old_status, new_status = old['status'].lower(), new['status'].lower()
    transitions['%s->%s' % (old_status, new_status)] += 1
    if old_status == 'concordant' and new_status == 'discordant':
        corrupted.append({'read': name, 'baseline_ps': old['PS'], 'run_ps': new['PS']})
    if old_status == 'discordant' and new_status == 'concordant':
        recovered.append({'read': name, 'baseline_ps': old['PS'], 'run_ps': new['PS']})
gained = [n for n in after if n not in before]
lost_by_ps = collections.Counter(before[n]['PS'] for n in lost)
lost_status = collections.Counter(before[n]['status'].lower() for n in lost)

summaries = {}
for tag, root in (('baseline', a.baseline), ('run', a.run)):
    summaries[tag] = json.load((root / 'eval/summary.json').open())

gaps_before, gaps_after = gap_status(tiers(a.baseline)), gap_status(tiers(a.run))
gap_moves = collections.Counter()
newly_joined, unjoined = [], []
for key, old in gaps_before.items():
    new = gaps_after.get(key, 'absent')
    if old != new:
        gap_moves['%s->%s' % (old, new)] += 1
        if new == 'joined' and old != 'joined':
            newly_joined.append(key)
        if old == 'joined' and new != 'joined':
            unjoined.append(key)

metrics = ('total_phased_reads', 'unphased_reads', 'total_phase_sets',
           'total_reads_evaluated', 'concordant_reads', 'discordant_reads',
           'hamming_error_rate', 'switch_errors', 'flip_errors',
           'phase_block_n50_bp', 'phase_block_aun_bp')
report = {
    'baseline': str(a.baseline),
    'run': str(a.run),
    'metrics': {m: {'baseline': summaries['baseline'][m], 'run': summaries['run'][m],
                    'delta': summaries['run'][m] - summaries['baseline'][m]} for m in metrics},
    'gate_concordant_to_discordant': len(corrupted),
    'discordant_to_concordant': len(recovered),
    'reads_gained': len(gained),
    'reads_lost': len(lost),
    'read_transitions': dict(transitions),
    'gap_status_changes': dict(gap_moves),
    'gaps_newly_joined': [list(k) for k in newly_joined],
    'gaps_no_longer_joined': [list(k) for k in unjoined],
    'gap_status_totals': {
        'baseline': dict(collections.Counter(gaps_before.values())),
        'run': dict(collections.Counter(gaps_after.values())),
    },
    'reads_lost_by_baseline_ps': dict(lost_by_ps.most_common(10)),
    'reads_lost_by_baseline_status': dict(lost_status),
    'corrupted_reads': corrupted[:50],
}
a.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
print(json.dumps({k: v for k, v in report.items()
                  if k not in ('corrupted_reads', 'read_transitions', 'gaps_newly_joined',
                               'gaps_no_longer_joined')}, indent=2, sort_keys=True))
