#!/usr/bin/env python3
"""Preserve correct reads and VCF evidence; explicitly audit the false graph SNP."""
from collections import Counter, defaultdict
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('reads', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
reads = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reads)
FALSE_SNP = ('CHM13#0#chr20', '25855631', '.', 'T', 'C')


def compare(before, after, truth):
    old, old_status = reads.assignments(before/'phased.bam', truth)
    new, new_status = reads.assignments(after/'phased.bam', truth)
    assert old.keys() <= new.keys(), (before, old.keys()-new.keys())
    lost = [q for q, state in old_status.items()
            if state == 'correct' and new_status.get(q) != 'correct']
    assert not lost, (before, lost)
    dropped = [q for q, (hp, ps) in old.items() if hp in (1, 2) and ps > 0 and
               (new[q][0] not in (1, 2) or new[q][1] <= 0)]
    assert not dropped, (before, dropped)
    old_counts, new_counts = Counter(old_status.values()), Counter(new_status.values())
    assert new_counts['discordant'] <= old_counts['discordant'], (before, new_counts)
    # Incorrect old assignments may change individually. Correct old molecules
    # must retain one coherent gauge for each old/new phase-set pair.
    gauges = defaultdict(set)
    for q, state in old_status.items():
        if state == 'correct':
            hp, ps = old[q]
            nhp, nps = new[q]
            gauges[(ps, nps)].add(hp != nhp)
    assert all(len(flips) == 1 for flips in gauges.values()), (before, gauges)
    old_vcf = reads.variants(before/'phased.vcf')
    new_vcf = reads.variants(after/'phased.vcf')
    removed = old_vcf.keys()-new_vcf.keys()
    assert removed <= {FALSE_SNP}, (before, removed)
    added = new_vcf.keys()-old_vcf.keys()
    assert all(25844279 <= int(key[1]) <= 25855615 for key in added), (before, added)
    changed_vcf = 0
    for key in old_vcf.keys() & new_vcf.keys():
        info, sample = old_vcf[key]
        ninfo, nsample = new_vcf[key]
        assert info == ninfo, (before, key)
        assert {k: v for k, v in sample.items() if k not in ('GT', 'PS')} == {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}, (before, key)
        gt, ngt = sample.get('GT'), nsample.get('GT')
        if gt != ngt:
            assert '|' in gt and ngt == '|'.join(reversed(gt.split('|'))), (before, key, gt, ngt)
            assert gauges.get((int(sample['PS']), int(nsample['PS']))) == {True}, (before, key)
        changed_vcf += sample != nsample
    return {'old_primary_output_reads': len(old), 'new_primary_output_reads': len(new),
            'old_counts': dict(old_counts), 'new_counts': dict(new_counts),
            'lost_correct_assignments': len(lost), 'lost_phased_assignments': len(dropped),
            'old_variant_records': len(old_vcf), 'new_variant_records': len(new_vcf),
            'removed_physically_contradicted_snp': list(removed),
            'added_recovery_variants': len(added),
            'added_variant_extent': [min(int(k[1]) for k in added), max(int(k[1]) for k in added)] if added else None,
            'retained_variant_alleles_counts_and_filters_preserved': True,
            'changed_variant_gauges': changed_vcf,
            'changed_read_tags': sum(old[q] != new[q] for q in old)}
