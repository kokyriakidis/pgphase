#!/usr/bin/env python3
"""Preserve prior native replay assignments and variants across the gap fix."""
import argparse
from collections import Counter, defaultdict
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--after', type=Path, default=Path('/tmp/pgphase-window-tests'))
args = parser.parse_args()
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}


def compare(before, after):
    old, old_status = audit.assignments(before/'phased.bam', truth)
    new, new_status = audit.assignments(after/'phased.bam', truth)
    assert old.keys() <= new.keys(), old.keys() - new.keys()
    assert all(new_status.get(q) == 'correct' for q, status in old_status.items()
               if status == 'correct'), before
    dropped = [q for q, (hp, ps) in old.items()
               if hp in (1, 2) and ps > 0 and (new[q][0] not in (1, 2) or new[q][1] <= 0)]
    # Withdrawing an unsupported discordant tag improves accuracy; report it
    # explicitly rather than treating it as preservation of a phased read.
    assert all(old_status.get(q) == 'discordant' and new_status.get(q) == 'unphased'
               for q in dropped), (before, dropped)
    old_counts, new_counts = Counter(old_status.values()), Counter(new_status.values())
    assert new_counts['discordant'] <= old_counts['discordant'], before
    changed = [q for q in old if old[q] != new[q]]
    gauge_maps = defaultdict(set)
    for q in changed:
        hp, ps = old[q]
        nhp, nps = new[q]
        if hp in (1, 2) and ps > 0 and nhp in (1, 2) and nps > 0:
            gauge_maps[(ps, nps)].add(hp != nhp)
    assert all(len(flips) == 1 for flips in gauge_maps.values()), (before, gauge_maps)
    old_vcf, new_vcf = audit.variants(before/'phased.vcf'), audit.variants(after/'phased.vcf')
    assert old_vcf.keys() <= new_vcf.keys(), before
    changed_vcf = 0
    for key, (info, sample) in old_vcf.items():
        ninfo, nsample = new_vcf[key]
        assert info == ninfo, (before, key)
        assert {k: v for k, v in sample.items() if k not in ('GT', 'PS')} == {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}, (before, key)
        gt, ngt = sample.get('GT'), nsample.get('GT')
        if gt != ngt:
            assert '|' in gt and ngt == '|'.join(reversed(gt.split('|'))), (before, key)
            assert gauge_maps.get((int(sample['PS']), int(nsample['PS']))) == {True}, (before, key)
        changed_vcf += sample != nsample
    return {'old_primary_output_reads': len(old), 'new_primary_output_reads': len(new),
            'old_correct': old_counts['correct'], 'new_correct': new_counts['correct'],
            'old_discordant': old_counts['discordant'], 'new_discordant': new_counts['discordant'],
            'old_variant_records': len(old_vcf), 'new_variant_records': len(new_vcf),
            'changed_variant_gauges': changed_vcf, 'changed_read_tags': len(changed),
            'lost_correct_assignments': 0, 'lost_phased_assignments': len(dropped),
            'discordant_tags_now_abstaining': dropped,
            'variant_alleles_counts_and_filters_preserved': True,
            'new_assignments': [q for q, (hp, ps) in new.items() if hp in (1, 2) and ps > 0 and
                                (q not in old or old[q][0] not in (1, 2))]}


reports = []
missing = []
for before_bam in sorted(Path('/tmp/pgphase71-final').rglob('phased.bam')):
    relative = before_bam.relative_to(Path('/tmp/pgphase71-final'))
    after_bam = args.after/relative
    if not after_bam.exists():
        missing.append(str(relative))
        continue
    result = compare(before_bam.parent, after_bam.parent)
    reports.append({'baseline': str(before_bam), 'output': str(relative), **result})
assert len(reports) == 219, f'Expected 219 prior native replays, got {len(reports)}'
assert not missing, missing
result = {'native_replays_compared': len(reports), 'unmatched_prior_replays': missing,
          'variant_alleles_counts_and_filters_preserved': True,
          'correct_assignments_preserved': True,
          'phased_assignments_preserved': not any(r['lost_phased_assignments'] for r in reports),
          'discordant_tags_now_abstaining': [
              {'output': r['output'], 'qname': q}
              for r in reports for q in r['discordant_tags_now_abstaining']],
          'changed_outputs': [r for r in reports if r['changed_read_tags'] or r['new_assignments']]}
(OUT/'panel-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_outputs'}, indent=2))
