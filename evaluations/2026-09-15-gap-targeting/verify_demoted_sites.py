#!/usr/bin/env python3
"""Ask whether the graph sites pgphase discards inside gaps carry real phase signal.

`add_demoted_site_linkage.py` showed that repeat-demoted het indels would link most
unresolved gaps. Linkage from a site is worth nothing unless the site's per-read
alleles actually follow the haplotypes, so each candidate is genotyped from the
alignment and scored against the read-level truth: a site whose alleles segregate
with truth is real signal that screening removed, one that partitions reads at
random is the noise the screen exists to catch.

Repeat-context indels are placed arbitrarily within their repeat run by the
aligner, so the allele is read from the net length change over a window rather
than from an exact anchor position -- an anchor test scores clean het indels at
2% and is useless here. The same scoring is applied to the sites the pipeline
trusts (CLEAN_HET_INDEL) and the ones its own MSA step verifies
(NOISY_CAND_HET), giving both a control and the precision of the verifier that
would gate the demoted sites. Finally, each gap's bottleneck linkage is
recomputed using only truth-validated sites, which is the honest size of the
opportunity.
"""
import argparse
import bisect
import collections
import csv
import random
import statistics
from pathlib import Path

import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--targets', type=Path, required=True)
p.add_argument('--graph-variants', type=Path, required=True)
p.add_argument('--hybrid-candidates', type=Path, required=True)
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--min-mapq', type=int, default=30)
p.add_argument('--min-truth-reads', type=int, default=10)
p.add_argument('--min-allele-reads', type=int, default=3)
p.add_argument('--informative-at', type=float, default=0.90)
p.add_argument('--phantom-below', type=float, default=0.70)
p.add_argument('--min-link-reads', type=int, default=2)
p.add_argument('--control-sample', type=int, default=200)
p.add_argument('--seed', type=int, default=5)
p.add_argument('--site-output', type=Path, required=True)
p.add_argument('--gap-output', type=Path, required=True)
a = p.parse_args()

BAM_CINS, BAM_CDEL, BAM_CREF_SKIP = 1, 2, 3
CONSUME_REF = frozenset((0, 2, 3, 7, 8))

truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap
bam = pysam.AlignmentFile(str(a.bam), 'rb')


def net_indel(read, lo, hi):
    delta, ref = 0, read.reference_start
    for op, length in read.cigartuples or ():
        if op == BAM_CINS:
            if lo <= ref <= hi:
                delta += length
        elif op in (BAM_CDEL, BAM_CREF_SKIP):
            if ref + length >= lo and ref <= hi:
                delta -= length
            ref += length
        elif op in CONSUME_REF:
            ref += length
    return delta


def genotype(pos, ref_allele, alt_allele, vtype):
    """read name -> 'ref'|'alt' at one site, read from the alignment."""
    calls = {}
    if vtype == 'SNP':
        for column in bam.pileup(a.contig, pos - 1, pos, truncate=True, min_base_quality=0,
                                 max_depth=10000, stepper='nofilter'):
            for read in column.pileups:
                if read.is_del or read.query_position is None:
                    continue
                base = read.alignment.query_sequence[read.query_position].upper()
                if base == alt_allele[:1].upper():
                    calls[read.alignment.query_name] = 'alt'
                elif base == ref_allele[:1].upper():
                    calls[read.alignment.query_name] = 'ref'
        return calls
    # TSV allele convention (src/collect_output.cpp): a deletion carries the
    # deleted bases in REF with ALT '.', an insertion carries the inserted
    # sequence in ALT against a single anchor base in REF. Comparing raw string
    # lengths silently skips every deletion.
    first_alt = alt_allele.split(',')[0]
    if vtype == 'DEL' or first_alt == '.':
        expect = -len(ref_allele)
    elif vtype == 'INS':
        expect = len(first_alt) - (len(ref_allele) - 1 if len(ref_allele) > 1 else 0)
    else:
        expect = len(first_alt) - len(ref_allele)
    if expect == 0:
        return {}
    span = max(len(ref_allele), 1)
    width = max(20, 2 * abs(expect) + 10)
    lo, hi = pos - width, pos + span + width
    for read in bam.fetch(a.contig, max(0, pos - 1), pos + span):
        if read.is_unmapped or read.cigartuples is None:
            continue
        if read.reference_start > pos - 1 or read.reference_end < pos + span:
            continue
        delta = net_indel(read, lo, hi)
        calls[read.query_name] = 'alt' if abs(delta - expect) < abs(delta) else 'ref'
    return calls


def score(calls):
    pairs = [(v, truth[k]) for k, v in calls.items() if k in truth]
    if len(pairs) < a.min_truth_reads:
        return None
    counts = collections.Counter(pairs)
    n_alt = sum(n for (allele, _), n in counts.items() if allele == 'alt')
    if min(n_alt, len(pairs) - n_alt) < a.min_allele_reads:
        return None
    straight = counts[('ref', 'MATERNAL')] + counts[('alt', 'PATERNAL')]
    flipped = counts[('ref', 'PATERNAL')] + counts[('alt', 'MATERNAL')]
    return {'reads': len(pairs), 'alt_fraction': n_alt / len(pairs),
            'segregation': max(straight, flipped) / len(pairs)}


by_category = collections.defaultdict(list)
with a.graph_variants.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        by_category[r['CATEGORY']].append((int(r['POS']), r['REF'], r['ALT'], r['TYPE']))
hybrid = collections.defaultdict(list)
with a.hybrid_candidates.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        hybrid[r['CATEGORY']].append((int(r['POS']), r['REF'], r['ALT'], r['TYPE']))
for table in (by_category, hybrid):
    for key in table:
        table[key].sort()

clean_positions = sorted(p for p, *_ in by_category['CLEAN_HET_SNP'] + by_category['CLEAN_HET_INDEL'])
targets = [r for r in csv.DictReader(a.targets.open(), delimiter='\t')
           if r['revised_cause'] == 'repeat_indels_would_bridge']


def in_window(sites, lo, hi):
    positions = [s[0] for s in sites]
    return sites[bisect.bisect_left(positions, lo):bisect.bisect_right(positions, hi)]


def linking_reads(reads, sites, cut):
    n = 0
    for start, end, mapq in reads:
        if mapq < a.min_mapq or start >= cut or end <= cut:
            continue
        left = bisect.bisect_right(sites, cut - 1) - bisect.bisect_left(sites, start)
        right = bisect.bisect_right(sites, end) - bisect.bisect_left(sites, cut + 1)
        if left and right:
            n += 1
    return n


site_rows, gap_rows = [], []
for gap in targets:
    lo, hi = int(gap['target_beg']), int(gap['target_end'])
    cut = int(gap['bottleneck_pos'])
    verdicts = collections.Counter()
    informative = []
    for pos, ref_allele, alt_allele, vtype in in_window(by_category['REP_HET_INDEL'], lo, hi):
        s = score(genotype(pos, ref_allele, alt_allele, vtype))
        if s is None:
            verdicts['unscorable'] += 1
            continue
        verdict = ('informative' if s['segregation'] >= a.informative_at
                   else 'phantom' if s['segregation'] < a.phantom_below else 'ambiguous')
        verdicts[verdict] += 1
        if verdict == 'informative':
            informative.append(pos)
        site_rows.append({'site_class': 'demoted_repeat_indel', 'pos': pos,
                          'gap_left': gap['gap_left'], 'gap_right': gap['gap_right'],
                          'verdict': verdict, **s})
    reads = [(r.reference_start, r.reference_end, r.mapping_quality)
             for r in bam.fetch(a.contig, max(0, lo), hi)
             if not r.is_unmapped and r.reference_end is not None]
    clean_here = clean_positions[bisect.bisect_left(clean_positions, lo):
                                 bisect.bisect_right(clean_positions, hi)]
    gap_rows.append({
        'gap_left': gap['gap_left'], 'gap_right': gap['gap_right'], 'gap_bp': gap['gap_bp'],
        'status': gap['status'], 'bottleneck_pos': cut,
        'demoted_informative': verdicts['informative'], 'demoted_ambiguous': verdicts['ambiguous'],
        'demoted_phantom': verdicts['phantom'], 'demoted_unscorable': verdicts['unscorable'],
        'link_clean': linking_reads(reads, clean_here, cut),
        'link_all_demoted': int(gap['link_with_repeat_indels']),
        'link_truth_validated': linking_reads(reads, sorted(clean_here + informative), cut),
    })

random.seed(a.seed)
for label, pool in (('control_clean_het_indel', hybrid['CLEAN_HET_INDEL']),
                    ('control_msa_verified_het', hybrid['NOISY_CAND_HET'])):
    for pos, ref_allele, alt_allele, vtype in random.sample(pool, min(a.control_sample, len(pool))):
        s = score(genotype(pos, ref_allele, alt_allele, vtype))
        if s is None:
            continue
        verdict = ('informative' if s['segregation'] >= a.informative_at
                   else 'phantom' if s['segregation'] < a.phantom_below else 'ambiguous')
        site_rows.append({'site_class': label, 'pos': pos, 'gap_left': '', 'gap_right': '',
                          'verdict': verdict, **s})

for path, rows, fields in (
        (a.site_output, site_rows, ['site_class', 'pos', 'gap_left', 'gap_right', 'verdict',
                                    'reads', 'alt_fraction', 'segregation']),
        (a.gap_output, gap_rows, list(gap_rows[0]))):
    with path.open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            r = dict(r)
            for k in ('alt_fraction', 'segregation'):
                if k in r:
                    r[k] = round(r[k], 4)
            w.writerow(r)

print('site verdicts by class:')
for cls in sorted({r['site_class'] for r in site_rows}):
    sel = [r for r in site_rows if r['site_class'] == cls]
    counts = collections.Counter(r['verdict'] for r in sel)
    print('  %-26s n=%-4d median seg %.3f   informative %3.0f%%  ambiguous %3.0f%%  phantom %3.0f%%' % (
        cls, len(sel), statistics.median([r['segregation'] for r in sel]),
        100 * counts['informative'] / len(sel), 100 * counts['ambiguous'] / len(sel),
        100 * counts['phantom'] / len(sel)))
kept = [r for r in gap_rows if r['link_truth_validated'] >= a.min_link_reads]
print()
print('gaps re-tested: %d' % len(gap_rows))
print('  linkage from all demoted indels        : %d gaps' % sum(
    1 for r in gap_rows if r['link_all_demoted'] >= a.min_link_reads))
print('  linkage from truth-VALIDATED sites only: %d gaps, %.2f Mb' % (
    len(kept), sum(int(r['gap_bp']) for r in kept) / 1e6))
print('  demoted sites in windows: %d informative, %d ambiguous, %d phantom, %d unscorable' % (
    sum(r['demoted_informative'] for r in gap_rows), sum(r['demoted_ambiguous'] for r in gap_rows),
    sum(r['demoted_phantom'] for r in gap_rows), sum(r['demoted_unscorable'] for r in gap_rows)))
