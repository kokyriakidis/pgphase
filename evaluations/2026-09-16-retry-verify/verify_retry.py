#!/usr/bin/env python3
"""Verify one window's phasing so that no check can pass by being blind.

Every wrong conclusion this file exists to prevent came from a check that
skipped something silently instead of failing:

  * A site was read as "phased and usable" from the candidate TSV's PHASE_SET
    column while the VCF emitted it `1|1` -- a candidate carries a phase set
    internally even when the solve put the alt allele on both haplotypes, and
    such a site links nothing. CATEGORY and emitted GENOTYPE are therefore
    reported as separate columns, and a category-het/genotype-hom site is an
    explicit finding.
  * A spacing was called bridgeable because sites existed on both sides,
    without counting the reads that cover both. Every consecutive pair of
    emitted het sites is counted here, and a link across a pair with no
    spanning read FAILS regardless of what the read gate says.
  * The read-level gate reported zero flipped reads for a block that was
    switched, because no read spans the switch point. A gate result is
    therefore never a verdict on its own: `gate_blind` is set whenever an
    unsupported link or an unscorable site exists in the window, and a PASS
    requires it to be false.
  * A per-site truth check silently dropped the one site that mattered because
    it fell below a read threshold. Unscorable sites are counted, listed and
    fail the verdict rather than being omitted.

Usage:
    verify_retry.py --vcf run/native.vcf --bam surjected.bam \
        --truth-map truth_hap.tsv --region-left L --region-right R \
        [--phased-bam run/phased.bam] [--baseline-phased-bam base/phased.bam] \
        [--competitor-vcf hiphase.vcf.gz --competitor-contig chr20] \
        --output verdict.json
"""
import argparse
import collections
import json
import re
import subprocess
from pathlib import Path

CONTIG = 'CHM13#0#chr20'
CONSUME_REF = frozenset('MDN=X')


def run(cmd):
    """Run a command and fail loudly. A silent non-zero exit has produced a
    wrong conclusion here before (region queries against an unaligned BAM)."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f'{" ".join(cmd)} exited {p.returncode}: {p.stderr[:400]}')
    return p.stdout


def ref_end(start, cigar):
    end = start
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in CONSUME_REF:
            end += int(n)
    return end


def load_vcf(path, lo, hi):
    """Every record in the window, classified. Nothing is dropped silently."""
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            pos = int(c[1])
            if not (lo <= pos <= hi):
                continue
            keys = c[8].split(':')
            vals = c[9].split(':')
            gt = vals[0]
            ps = vals[keys.index('PS')].strip() if 'PS' in keys else None
            cat = ''
            for field in c[7].split(';'):
                if field.startswith('CAT='):
                    cat = field[4:]
            phased = '|' in gt
            hom = gt in ('0|0', '1|1', '0/0', '1/1')
            out.append(dict(pos=pos, ref=c[3], alt=c[4], gt=gt, ps=ps, cat=cat,
                            phased=phased, hom=hom,
                            usable=phased and not hom and ps not in (None, '', '0', '.')))
    return out


def load_reads(lo, hi, bam, min_mapq=1):
    reads = {}
    for line in run(['samtools', 'view', '-q', str(min_mapq), bam,
                     f'{CONTIG}:{lo}-{hi}']).splitlines():
        f = line.split('\t')
        if int(f[1]) & 0x900:
            continue
        reads[f[0]] = (int(f[3]), f[5], f[9])
    return reads


def spanning(reads, a, b):
    return sum(1 for st, cig, _ in reads.values() if st <= a and ref_end(st, cig) >= b)


def read_allele(read, pos, ref, alt, flank=25):
    """Substitution: the base at the position. Indel: net insert-minus-delete
    length over a window, never the anchor position -- the aligner places a
    repeat-context indel arbitrarily within its run."""
    start, cigar, seq = read
    ops = [(int(n), o) for n, o in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
    if len(ref) == 1 and len(alt) == 1:
        r = start
        q = 0
        for n, op in ops:
            if op in 'M=X':
                if r <= pos < r + n:
                    base = seq[q + pos - r]
                    return 0 if base == ref else (1 if base == alt else None)
                r += n
                q += n
            elif op in 'DN':
                r += n
            elif op in 'IS':
                q += n
        return None
    expect = len(alt) - len(ref)
    lo_, hi_ = pos - flank, pos + flank
    r = start
    delta = 0
    for n, op in ops:
        if op == 'I':
            if lo_ <= r <= hi_:
                delta += n
        elif op in 'DN':
            if r + n >= lo_ and r <= hi_:
                delta -= n
            r += n
        elif op in 'M=X':
            r += n
    if start > lo_ or r < hi_:
        return None
    return 1 if abs(delta - expect) < abs(delta) else 0


def site_orientation(reads, truth, rec, min_scored):
    """Which parent this site's first GT allele carries, or an explicit
    UNSCORABLE with the reason -- never a silent skip."""
    table = collections.Counter()
    for qname, read in reads.items():
        parent = truth.get(qname)
        if parent is None:
            continue
        allele = read_allele(read, rec['pos'], rec['ref'], rec['alt'])
        if allele is None:
            continue
        hap = 1 if rec['gt'].split('|')[0] == str(allele) else 2
        table[(hap, parent)] += 1
    n = sum(table.values())
    if n < min_scored:
        return dict(orientation=None, confidence=None, scored=n,
                    unscorable='too_few_scored_reads')
    mat = table[(1, 'MATERNAL')] + table[(2, 'PATERNAL')]
    pat = table[(1, 'PATERNAL')] + table[(2, 'MATERNAL')]
    if mat == pat:
        return dict(orientation=None, confidence=0.5, scored=n,
                    unscorable='tied_orientation')
    return dict(orientation='MAT' if mat > pat else 'PAT',
                confidence=max(mat, pat) / n, scored=n, unscorable=None)


def tag_map(bam):
    out = {}
    for line in run(['samtools', 'view', bam]).splitlines():
        f = line.split('\t')
        hap = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hap = int(x[5:])
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hap and ps:
            out[f[0]] = (hap, ps)
    return out


def concordance(tags, truth):
    by_ps = collections.defaultdict(list)
    for qname, (hap, ps) in tags.items():
        if qname in truth:
            by_ps[ps].append((qname, hap))
    out = {}
    for ps, members in by_ps.items():
        table = collections.Counter((hap, truth[q]) for q, hap in members)
        straight = table[(1, 'MATERNAL')] + table[(2, 'PATERNAL')]
        flipped = table[(1, 'PATERNAL')] + table[(2, 'MATERNAL')]
        flip = flipped > straight
        for q, hap in members:
            want = 'MATERNAL' if (hap == 1) != flip else 'PATERNAL'
            out[q] = truth[q] == want
    return out


p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--vcf', type=Path, required=True)
p.add_argument('--bam', type=Path, required=True, help='surjected input BAM')
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--region-left', type=int, required=True)
p.add_argument('--region-right', type=int, required=True)
p.add_argument('--context', type=int, default=60000,
               help='bp each side included so the flanking blocks are seen')
p.add_argument('--phased-bam', type=Path)
p.add_argument('--baseline-phased-bam', type=Path)
p.add_argument('--competitor-vcf', type=Path)
p.add_argument('--competitor-contig', default='chr20')
p.add_argument('--min-scored-reads', type=int, default=6,
               help='reads a site needs before its orientation is trusted')
p.add_argument('--switch-run', type=int, default=2,
               help='sites that must agree on each side before a switch is declared')
p.add_argument('--min-confidence', type=float, default=0.8,
               help='orientation confidence below this counts as unscorable')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

LO = max(1, a.region_left - a.context)
HI = a.region_right + a.context
truth = dict(l.rstrip('\n').split('\t') for l in open(a.truth_map))
reads = load_reads(LO, HI, str(a.bam))
records = load_vcf(a.vcf, LO, HI)

findings = []

# 1. Category says het, genotype says hom. Such a site is in the table, carries a
#    phase set, and links nothing.
mislabelled = [r for r in records
               if 'HET' in r['cat'] and r['hom']]
for r in mislabelled:
    findings.append(dict(kind='category_het_genotype_hom', pos=r['pos'],
                         cat=r['cat'], gt=r['gt'],
                         inside_region=a.region_left <= r['pos'] <= a.region_right))

# 2. The chain: consecutive usable het sites per phase set, with spanning reads.
usable = [r for r in records if r['usable']]
by_ps = collections.defaultdict(list)
for r in usable:
    by_ps[r['ps']].append(r)
links = []
for ps, members in by_ps.items():
    members.sort(key=lambda r: r['pos'])
    for i in range(len(members) - 1):
        left, right = members[i]['pos'], members[i + 1]['pos']
        n = spanning(reads, left, right)
        row = dict(ps=ps, left=left, right=right, bp=right - left, spanning_reads=n,
                   supported=n > 0,
                   inside_region=not (right < a.region_left or left > a.region_right))
        links.append(row)
        if n == 0:
            findings.append(dict(kind='unsupported_link', **row))

# 3. Per-site truth orientation, with unscorable sites named rather than dropped.
scored = []
for r in usable:
    o = site_orientation(reads, truth, r, a.min_scored_reads)
    if o['unscorable'] is None and o['confidence'] < a.min_confidence:
        o = dict(o, unscorable='below_min_confidence')
    scored.append(dict(pos=r['pos'], ps=r['ps'], cat=r['cat'], **o))
unscorable = [s for s in scored if s['unscorable']]
for s in unscorable:
    findings.append(dict(kind='unscorable_site', pos=s['pos'], ps=s['ps'],
                         reason=s['unscorable'], scored=s['scored'],
                         inside_region=a.region_left <= s['pos'] <= a.region_right))

# 4. Switch localisation: a sign change that PERSISTS, inside one phase set.
#
# Two things make a single-pair comparison unusable. Several records can sit at
# one position (multiallelic decomposition, or an indel written at adjacent
# anchors), and comparing them to each other reports a switch between a site and
# itself -- three of the six this check first emitted were of that kind. And one
# site whose orientation is called wrong produces two spurious sign changes. So
# positions are collapsed to a majority orientation first, and a switch is
# declared only where a run of at least --switch-run sites on each side disagrees.
switches = []
for ps, members in by_ps.items():
    by_pos = collections.defaultdict(collections.Counter)
    for s in scored:
        if s['ps'] == ps and s['orientation']:
            by_pos[s['pos']][s['orientation']] += 1
    seq = []
    for pos in sorted(by_pos):
        counts = by_pos[pos]
        if len(counts) > 1 and counts.most_common(1)[0][1] == list(counts.values())[1]:
            continue  # position disagrees with itself: no usable orientation
        seq.append((pos, counts.most_common(1)[0][0]))
    run_len = a.switch_run
    for i in range(len(seq) - 1):
        left_run = [o for _, o in seq[max(0, i - run_len + 1):i + 1]]
        right_run = [o for _, o in seq[i + 1:i + 1 + run_len]]
        if len(left_run) < run_len or len(right_run) < run_len:
            continue
        if len(set(left_run)) != 1 or len(set(right_run)) != 1:
            continue
        if left_run[0] == right_run[0]:
            continue
        gap_reads = spanning(reads, seq[i][0], seq[i + 1][0])
        row = dict(ps=ps, left=seq[i][0], right=seq[i + 1][0],
                   left_orientation=left_run[0], right_orientation=right_run[0],
                   spanning_reads=gap_reads, run=run_len)
        switches.append(row)
        findings.append(dict(kind='switch_within_phase_set', **row))

# 5. What a competitor uses here, and what we do at those positions.
competitor = None
if a.competitor_vcf:
    comp_rows = []
    for line in run(['tabix', str(a.competitor_vcf),
                     f'{a.competitor_contig}:{LO}-{HI}']).splitlines():
        c = line.split('\t')
        gt = c[9].split(':')[0]
        if '|' not in gt or gt in ('0|0', '1|1'):
            continue
        comp_rows.append((int(c[1]), c[3], c[4], gt))
    ours_by_pos = {r['pos']: r for r in records}
    census = collections.Counter()
    detail = []
    for pos, ref, alt, gt in comp_rows:
        hit = None
        for d in (0, 1, -1, 2, -2):
            if pos + d in ours_by_pos:
                hit = ours_by_pos[pos + d]
                break
        if hit is None:
            state = 'absent_from_our_vcf'
        elif hit['usable']:
            state = 'we_phase_it_het'
        elif hit['hom']:
            state = 'we_call_it_hom'
        else:
            state = 'we_leave_it_unphased'
        census[state] += 1
        if state != 'we_phase_it_het':
            detail.append(dict(pos=pos, state=state,
                               our_cat=hit['cat'] if hit else None,
                               our_gt=hit['gt'] if hit else None,
                               inside_region=a.region_left <= pos <= a.region_right))
            findings.append(dict(kind='competitor_site_unused', pos=pos, state=state,
                                 our_cat=hit['cat'] if hit else None,
                                 our_gt=hit['gt'] if hit else None,
                                 inside_region=a.region_left <= pos <= a.region_right))
    competitor = dict(phased_hets=len(comp_rows), census=dict(census), detail=detail)

# 6. The read-level gate -- reported, never trusted alone.
gate = None
if a.phased_bam and a.baseline_phased_bam:
    before = tag_map(str(a.baseline_phased_bam))
    after = tag_map(str(a.phased_bam))
    cb, ca = concordance(before, truth), concordance(after, truth)
    x = collections.Counter()
    for q in set(cb) | set(ca):
        x[(cb.get(q), ca.get(q))] += 1
    gate = dict(tagged_before=len(before), tagged_after=len(after),
                concordant_before=sum(cb.values()), concordant_after=sum(ca.values()),
                concordant_to_discordant=x[(True, False)],
                concordant_lost=x[(True, None)],
                newly_concordant=x[(None, True)], newly_discordant=x[(None, False)])

# 7. Verdict. A gate result cannot carry a PASS on its own.
in_region = lambda f: f.get('inside_region', True)
unsupported_in_region = [f for f in findings
                         if f['kind'] == 'unsupported_link' and in_region(f)]
unscorable_in_region = [f for f in findings
                        if f['kind'] == 'unscorable_site' and in_region(f)]
gate_blind = bool(unsupported_in_region or unscorable_in_region)
blocking = [f for f in findings if f['kind'] in
            ('unsupported_link', 'switch_within_phase_set', 'category_het_genotype_hom')
            and in_region(f)]
if gate and gate['concordant_to_discordant'] > 0:
    blocking.append(dict(kind='gate_flipped_reads',
                         n=gate['concordant_to_discordant']))
verdict = dict(region=[a.region_left, a.region_right],
               usable_het_sites=len(usable),
               usable_het_sites_in_region=sum(
                   1 for r in usable if a.region_left <= r['pos'] <= a.region_right),
               links=len(links), unsupported_links=len(unsupported_in_region),
               unscorable_sites=len(unscorable_in_region),
               switches=len(switches),
               category_het_genotype_hom=len([f for f in findings
                                              if f['kind'] == 'category_het_genotype_hom']),
               gate=gate, gate_blind=gate_blind,
               competitor=competitor,
               passes=not blocking and not gate_blind)

a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(dict(verdict=verdict, findings=findings,
                                    links=links, sites=scored), indent=1))

print('region %d-%d   usable het sites %d (%d inside)' % (
    a.region_left, a.region_right, verdict['usable_het_sites'],
    verdict['usable_het_sites_in_region']))
print('links %d, of which unsupported (0 spanning reads) %d' % (
    verdict['links'], verdict['unsupported_links']))
for f in [f for f in findings if f['kind'] == 'unsupported_link' and in_region(f)]:
    print('   UNSUPPORTED %d -> %d (%.1f kb, %d spanning reads)' % (
        f['left'], f['right'], f['bp'] / 1e3, f['spanning_reads']))
print('sites unscorable against truth: %d' % verdict['unscorable_sites'])
for f in unscorable_in_region:
    print('   UNSCORABLE %d (%s, %d scored reads)' % (f['pos'], f['reason'], f['scored']))
print('category-het / genotype-hom sites: %d' % verdict['category_het_genotype_hom'])
for f in [f for f in findings if f['kind'] == 'category_het_genotype_hom']:
    print('   MISLABELLED %d cat=%s gt=%s%s' % (
        f['pos'], f['cat'], f['gt'], '  (inside region)' if f['inside_region'] else ''))
print('switches within a phase set: %d' % verdict['switches'])
for s in switches:
    print('   SWITCH in PS=%s between %d (%s) and %d (%s), %d spanning reads' % (
        s['ps'], s['left'], s['left_orientation'], s['right'],
        s['right_orientation'], s['spanning_reads']))
if competitor:
    print('competitor phased hets in window: %d   %s' % (
        competitor['phased_hets'], competitor['census']))
    for d in competitor['detail']:
        if d['inside_region']:
            print('   UNUSED %d %s (our cat=%s gt=%s)' % (
                d['pos'], d['state'], d['our_cat'], d['our_gt']))
if gate:
    print('gate  tagged %d -> %d   concordant %d -> %d   conc->disc %d   lost %d   new %dc/%dd'
          % (gate['tagged_before'], gate['tagged_after'], gate['concordant_before'],
             gate['concordant_after'], gate['concordant_to_discordant'],
             gate['concordant_lost'], gate['newly_concordant'],
             gate['newly_discordant']))
    print('gate_blind: %s%s' % (gate_blind,
          '  (an unsupported link or unscorable site exists here, so the gate'
          ' cannot see a switch across it)' if gate_blind else ''))
print('VERDICT: %s' % ('PASS' if verdict['passes'] else 'FAIL'))
print('wrote %s' % a.output)
