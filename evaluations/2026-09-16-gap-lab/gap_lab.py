#!/usr/bin/env python3
"""One command per gap: gather the evidence, phase it, stitch it, score it.

Every gap diagnosis in this project has been assembled by hand, which is slow and
makes each experiment's verdict depend on which checks were remembered. This runs
the whole loop for one gap and prints a verdict:

  1. BASELINE -- the pipeline as shipped over the gap plus a flank, which is what
     any change has to beat. Its phase blocks define the gauge the gap's flanks
     are already in.
  2. EVIDENCE -- every candidate the BAM channel finds inside the gap interval,
     with its category, and (with --truth-map) how well each site's allele
     partition follows the read truth. A site that does not segregate is not
     evidence, whatever its category says.
  3. GAP PHASING -- the BAM channel run on the interval. This is the existing
     phasing machinery, not a reimplementation: whatever it produces for the gap
     is what a scoped change inside the pipeline would have to reproduce.
  4. STITCH -- link each gap block to the two flanks by read identity, using the
     project's own net-margin orientation standard, and compose the parity. A gap
     is CLOSED only when one composed frame reaches both flanks.
  5. GATE -- read-level scoring against truth before and after. A stitch that
     turns concordant reads discordant is a regression however many reads it
     adds, so the verdict is pass/fail on that, not on the join count.

Phase-set extents come from the phased VCF (first to last phased variant). Read
start positions understate a block's reach by up to a read length at each end and
must not be used for this.
"""
import argparse
import collections
import csv
import json
import subprocess
from pathlib import Path

REF = 'test_data/chm13v2.0.chr20.renamed.fa'
BAM = 'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'
SITES = 'test_data/chr20.sites.striped.vcf.gz'
GAF = 'test_data/HG002.chr20.annotated.coord.gaf.gz'
CONTIG = 'CHM13#0#chr20'


def run(cmd, log):
    """Run a pipeline arm, failing loudly. A silent non-zero exit here has
    produced wrong conclusions before: an empty output looks like 'no evidence'."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w') as stream:
        p = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise SystemExit(f'arm failed ({p.returncode}): {" ".join(cmd[:3])}... see {log}')


def vcf_phase_sets(path):
    """ps -> (first phased variant, last phased variant, site count)."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or '|' not in c[9].split(':')[0]:
                continue
            keys = c[8].split(':')
            if 'PS' not in keys:
                continue
            ps = c[9].split(':')[keys.index('PS')].strip()
            pos = int(c[1])
            lo, hi, n = out.get(ps, (pos, pos, 0))
            out[ps] = (min(lo, pos), max(hi, pos), n + 1)
    return out


def bam_tags(path, region=None):
    """qname -> (hap, ps). Only reads the aligner placed and the pipeline tagged.

    A chromosome-wide phased BAM is read by region when an index exists: the
    gauge that matters is the blocks around this gap, and the gate is scored over
    the same read population.
    """
    cmd = ['samtools', 'view', str(path)]
    if region and Path(str(path) + '.bai').exists():
        cmd.append(region)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f'samtools view failed on {path}')
    out = {}
    for line in p.stdout.splitlines():
        f = line.split('\t')
        hp = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hp = int(x[5:])
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hp and ps:
            out[f[0]] = (hp, ps)
    return out


def votes(a_tags, a_ps, b_tags, b_ps):
    """n11/n12/n21/n22 over reads carrying a_ps in one labelling and b_ps in the other."""
    t = [0, 0, 0, 0]
    shared = 0
    for q, (ha, pa) in a_tags.items():
        if pa != a_ps:
            continue
        hb_ps = b_tags.get(q)
        if hb_ps is None or hb_ps[1] != b_ps:
            continue
        shared += 1
        t[(ha - 1) * 2 + (hb_ps[0] - 1)] += 1
    return t, shared


def decide(t, margin):
    """The project's net-margin rule: merge iff |(n12+n21)-(n11+n22)| > margin."""
    straight = t[0] + t[3]
    flipped = t[1] + t[2]
    if abs(flipped - straight) <= margin:
        return None
    return flipped > straight


def score(tags, truth):
    """Per block, orient to whichever assignment truth prefers, then per-read correctness."""
    by_ps = collections.defaultdict(list)
    for q, (h, ps) in tags.items():
        if q in truth:
            by_ps[ps].append((q, h))
    correct = {}
    for ps, members in by_ps.items():
        c = collections.Counter((h, truth[q]) for q, h in members)
        good = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
        bad = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
        flip = bad > good
        for q, h in members:
            expect = 'MATERNAL' if (h == 1) != flip else 'PATERNAL'
            correct[q] = truth[q] == expect
    return correct


def segregation(pos, ref, alt, vtype, flank, truth):
    """Read a site's alleles from the alignment and score them against read truth.

    Indels are read from the net insert-minus-delete length over a window: the
    aligner places a repeat-context indel arbitrarily within its run, so the
    anchor position says nothing. Substitutions are read as the base at the site
    (SAM POS and VCF POS are both 1-based and the reference walk starts at POS).
    """
    span = max(len(ref), 1)
    lo, hi = pos - flank, pos + span + flank
    p = subprocess.run(['samtools', 'view', '-q', '1', BAM,
                        f'{CONTIG}:{max(1, lo - 50)}-{hi + 50}'],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return None, 0
    calls = {}
    for line in p.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 10:
            continue
        name, start, cig, seq = f[0], int(f[3]), f[5], f[9]
        ops, num = [], ''
        for ch in cig:
            if ch.isdigit():
                num += ch
            else:
                ops.append((int(num), ch))
                num = ''
        if vtype == 'SNP':
            refp, qi, base = start, 0, None
            for length, op in ops:
                if op in 'M=X':
                    if refp <= pos < refp + length:
                        base = seq[qi + (pos - refp)]
                        break
                    refp += length
                    qi += length
                elif op == 'I':
                    qi += length
                elif op in 'DN':
                    refp += length
                elif op == 'S':
                    qi += length
            if base is None:
                continue
            if base.upper() == alt.upper():
                calls[name] = 1
            elif base.upper() == ref.upper():
                calls[name] = 0
        else:
            first_alt = alt.split(',')[0]
            expect = (-len(ref) if first_alt == '.' else len(first_alt) - len(ref))
            if expect == 0:
                continue
            refp, delta, rs = start, 0, start
            for length, op in ops:
                if op == 'I':
                    if lo <= refp <= hi:
                        delta += length
                elif op in 'DN':
                    if refp + length >= lo and refp <= hi:
                        delta -= length
                    refp += length
                elif op in 'M=X':
                    refp += length
            if rs > lo or refp < hi:
                continue
            if abs(delta - expect) < abs(delta):
                calls[name] = 1
            elif abs(delta) < abs(delta - expect):
                calls[name] = 0
    scored = [(v, truth[q]) for q, v in calls.items() if q in truth]
    if len(scored) < 10:
        return None, len(scored)
    c = collections.Counter(scored)
    straight = c[(1, 'MATERNAL')] + c[(0, 'PATERNAL')]
    flipped = c[(1, 'PATERNAL')] + c[(0, 'MATERNAL')]
    return max(straight, flipped) / len(scored), len(scored)



def vcf_sites(path):
    """pos -> (ref, alt, gt, ps) for phased heterozygous records."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            gt = c[9].split(':')[0]
            if '|' not in gt or gt in ('0|0', '1|1'):
                continue
            keys = c[8].split(':')
            if 'PS' not in keys:
                continue
            out[int(c[1])] = (c[3], c[4], gt, c[9].split(':')[keys.index('PS')].strip())
    return out


def load_reads(lo, hi):
    """qname -> (start, cigar, seq) for reads overlapping a window."""
    p_ = subprocess.run(['samtools', 'view', '-q', '1', BAM, f'{CONTIG}:{lo}-{hi}'],
                        capture_output=True, text=True)
    if p_.returncode != 0:
        raise SystemExit('samtools view failed while loading bridge reads')
    out = {}
    for line in p_.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 10:
            continue
        out[f[0]] = (int(f[3]), f[5], f[9])
    return out


def cigar_ops(cig):
    num = ''
    for ch in cig:
        if ch.isdigit():
            num += ch
        else:
            yield int(num), ch
            num = ''


def read_allele(read, pos, ref, alt, flank):
    """0 (ref), 1 (alt) or None, using the same conventions as segregation()."""
    start, cig, seq = read
    ops = list(cigar_ops(cig))
    if len(ref) == 1 and len(alt) == 1:
        refp, qi = start, 0
        for length, op in ops:
            if op in 'M=X':
                if refp <= pos < refp + length:
                    b = seq[qi + (pos - refp)].upper()
                    return 1 if b == alt.upper() else (0 if b == ref.upper() else None)
                refp += length
                qi += length
            elif op == 'I':
                qi += length
            elif op in 'DN':
                refp += length
            elif op == 'S':
                qi += length
        return None
    expect = len(alt) - len(ref)
    if expect == 0:
        return None
    lo_, hi_ = pos - flank, pos + max(len(ref), 1) + flank
    refp, delta = start, 0
    for length, op in ops:
        if op == 'I':
            if lo_ <= refp <= hi_:
                delta += length
        elif op in 'DN':
            if refp + length >= lo_ and refp <= hi_:
                delta -= length
            refp += length
        elif op in 'M=X':
            refp += length
    if start > lo_ or refp < hi_:
        return None
    if abs(delta - expect) < abs(delta):
        return 1
    if abs(delta) < abs(delta - expect):
        return 0
    return None


def bridge_vote(left_sites, right_sites, reads, seam_window, n_sites, min_obs,
                site_flank, seam_span=None):
    """Orient two blocks using reads that cross the seam, from ALLELES not PS tags.

    A read is tagged with at most one phase set, so two blocks that split at the
    same position share no tagged reads and cannot be oriented against each other
    by tag identity -- which is the state every gap examined here is in. The
    alleles are still there: a read crossing the seam observes sites belonging to
    both blocks, and each block's own phased genotypes say which haplotype each
    allele belongs to. That is the vote this builds.
    """
    # Take every site of each adjacent block that a bridge read could actually
    # reach, not a fixed count nearest the seam. A HiFi read crossing the seam
    # observes as many sites as fall inside it -- capping at eight discarded most
    # of them where the flank is site-dense (the right flank here holds 916
    # sites, about one per 800 bp), which weakened each read's own call and
    # dropped reads that would otherwise have cleared min_obs. Sites beyond read
    # reach cannot contribute however many of them the block holds, so the bound
    # is the read length, not the block.
    left = sorted(left_sites)
    right = sorted(right_sites)
    if not left or not right:
        return None
    if seam_span:
        seam_guess = (left[-1] + right[0]) // 2
        left = [pos for pos in left if pos >= seam_guess - seam_span]
        right = [pos for pos in right if pos <= seam_guess + seam_span]
        if not left or not right:
            return None
    left = left[-n_sites:]
    right = right[:n_sites]
    # A bridge read has to cross the SEAM -- the point between the blocks -- and
    # observe enough sites on each side of it. Requiring it to span every
    # selected site instead is unsatisfiable: the sites reach tens of kb back
    # into each block while the seam here is 220 bp wide, and the first version
    # of this vote reported zero crossing reads for exactly that reason.
    seam = (max(left) + min(right)) // 2

    def hap_of(read, positions, table):
        tab = collections.Counter()
        for pos in positions:
            ref, alt, gt, _ = table[pos]
            al = read_allele(read, pos, ref, alt, site_flank)
            if al is None:
                continue
            first = gt.split('|')[0]
            tab[1 if (first == str(al)) else 2] += 1
        if sum(tab.values()) < min_obs or tab[1] == tab[2]:
            return None
        return 1 if tab[1] > tab[2] else 2

    def ref_end(read):
        end = read[0]
        for length, op in cigar_ops(read[1]):
            if op in 'M=XDN':
                end += length
        return end

    t = [0, 0, 0, 0]
    crossing = 0
    for q, read in reads.items():
        if read[0] > seam - seam_window or ref_end(read) < seam + seam_window:
            continue
        crossing += 1
        ha = hap_of(read, left, left_sites)
        hb = hap_of(read, right, right_sites)
        if ha is None or hb is None:
            continue
        t[(ha - 1) * 2 + (hb - 1)] += 1
    return dict(votes=t, crossing=crossing, voters=sum(t), seam=seam,
                n_left_sites=len(left), n_right_sites=len(right))



def site_orientation(positions, table, reads, truth, site_flank, min_obs=10):
    """Which truth haplotype does hap1 carry at these sites, per the genotypes?

    A flank with no tagged reads contributes nothing to the read-level gate, so a
    wrong orientation on that side is invisible to it: every read in the merged
    block comes from the other side and keeps its own relative labelling. The
    genotypes are still checkable -- each site's alleles can be read off the
    alignment and compared with the read truth -- so the LINK is validated here
    even when the flank's block cannot be.
    """
    votes_ = collections.Counter()
    used = 0
    for pos in positions:
        ref, alt, gt, _ = table[pos]
        tab = collections.Counter()
        for q, read in reads.items():
            if q not in truth:
                continue
            al = read_allele(read, pos, ref, alt, site_flank)
            if al is None:
                continue
            first = gt.split('|')[0]
            hap = 1 if (first == str(al)) else 2
            tab[(hap, truth[q])] += 1
        n = sum(tab.values())
        if n < min_obs:
            continue
        mat = tab[(1, 'MATERNAL')] + tab[(2, 'PATERNAL')]
        pat = tab[(1, 'PATERNAL')] + tab[(2, 'MATERNAL')]
        if mat == pat:
            continue
        used += 1
        votes_['MATERNAL' if mat > pat else 'PATERNAL'] += 1
    if not used:
        return None, 0, votes_
    winner = 'MATERNAL' if votes_['MATERNAL'] >= votes_['PATERNAL'] else 'PATERNAL'
    return winner, used, votes_


p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--gap-left', type=int, required=True)
p.add_argument('--gap-right', type=int, required=True)
p.add_argument('--flank', type=int, default=50000,
               help='window for the BASELINE arm when no --baseline-vcf is given')
p.add_argument('--baseline-vcf', type=Path,
               help='phased VCF from the real pipeline run. The flanks a gap is '
                    'stitched to are the frozen chromosome-wide blocks, not blocks '
                    'from a windowed re-run: re-running over gap+flank re-solves '
                    'the flanks and can move the very boundaries under test.')
p.add_argument('--baseline-bam', type=Path,
               help='phased BAM from the same run (read tags define the gauge)')
p.add_argument('--gauge-window', type=int, default=200000,
               help='bp each side of the gap read from the baseline BAM for the gauge')
p.add_argument('--gap-margin', type=int, default=5000,
               help='read context around the gap for the GAP arm. The machinery '
                    'must run on the gap interval so its blocks are gap-local; '
                    'running it over gap+50kb inherits the window\'s own breaks '
                    'and produces blocks that are not the gap\'s phasing at all.')
p.add_argument('--threads', type=int, default=8)
p.add_argument('--min-mapq', type=int, default=1)
p.add_argument('--stitch-margin', type=int, default=2)
p.add_argument('--seam-read-window', type=int, default=60000,
               help='bp each side of the gap from which bridge reads are loaded')
p.add_argument('--seam-window', type=int, default=2000,
               help='bp a bridge read must extend past the outermost seam site')
p.add_argument('--seam-sites', type=int, default=200,
               help='cap on sites per side for the allele-level bridge vote')
p.add_argument('--seam-span', type=int, default=30000,
               help='bp each side of the seam from which block sites are taken; '
                    'set to a read length, since sites a bridge read cannot '
                    'reach contribute nothing however many the block holds')
p.add_argument('--min-site-obs', type=int, default=2,
               help='site observations a bridge read needs on each side to vote')
p.add_argument('--max-new-discordant', type=int, default=0,
               help='newly tagged discordant reads tolerated before the gate fails')
p.add_argument('--min-flank-reads', type=int, default=5,
               help='tagged reads a baseline block needs before it counts as a flank')
p.add_argument('--site-flank', type=int, default=25)
p.add_argument('--truth-map', type=Path)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--out', type=Path, help='verdict JSON')
p.add_argument('--reuse', action='store_true', help='skip arms whose outputs exist')
a = p.parse_args()

GL, GR = a.gap_left, a.gap_right
lo, hi = GL - a.flank, GR + a.flank
region = f'{CONTIG}:{lo}-{hi}'
a.work.mkdir(parents=True, exist_ok=True)
truth = {}
if a.truth_map:
    truth = dict(line.rstrip('\n').split('\t') for line in a.truth_map.open())

# ---- 1. BASELINE: the frozen gauge. Prefer the real pipeline run's outputs;
# only fall back to a windowed re-run when none is supplied.
base = a.work / 'baseline'
base.mkdir(exist_ok=True)
if a.baseline_vcf and a.baseline_bam:
    base_vcf, base_bam = a.baseline_vcf, a.baseline_bam
    baseline_source = f'pipeline run ({base_vcf})'
else:
    base_vcf, base_bam = base / 'native.vcf', base / 'phased.bam'
    baseline_source = f'windowed re-run ({region})'
    if not (a.reuse and base_vcf.exists()):
        run(['./pgphase', 'collect-hybrid-variation', '--ref', REF, '--bam', BAM,
             '--graph-sites', SITES, '--gaf', GAF, '-r', region,
             '-t', str(a.threads), '-q', str(a.min_mapq),
             '--link-by-alleles', '--block-link-window', '8',
             '--min-read-margin', '2',
             '--recover-gaps', '--gap-recovery-report', str(base / 'tiers.tsv'),
             '-o', str(base / 'candidates.tsv'),
             '--phased-vcf-out', str(base_vcf),
             '-b', str(base_bam)], base / 'run.log')

# ---- 2/3. GAP PHASING: the machinery ON THE GAP INTERVAL, plus only enough
# margin for read context. Its blocks are then the gap's own phasing, which is
# what gets stitched -- the same shape as recovery's proposal.
gap_region = f'{CONTIG}:{max(1, GL - a.gap_margin)}-{GR + a.gap_margin}'
gapdir = a.work / 'gap_bam'
gapdir.mkdir(exist_ok=True)
if not (a.reuse and (gapdir / 'native.vcf').exists()):
    run(['./pgphase', 'collect-bam-variation', '--ref', REF, '--bam', BAM,
         '-r', gap_region, '-t', str(a.threads), '-q', str(a.min_mapq),
         '-o', str(gapdir / 'candidates.tsv'),
         '--phased-vcf-out', str(gapdir / 'native.vcf'),
         '-b', str(gapdir / 'phased.bam')], gapdir / 'run.log')

base_sets = vcf_phase_sets(base_vcf)
gap_sets = vcf_phase_sets(gapdir / 'native.vcf')
gauge_region = f'{CONTIG}:{max(1, GL - a.gauge_window)}-{GR + a.gauge_window}'
base_tags = bam_tags(base_bam, gauge_region)
gap_tags = bam_tags(gapdir / 'phased.bam')
print(f'  baseline gauge: {baseline_source}')
print(f'  gap arm region: {gap_region} (gap + {a.gap_margin} bp read context)')

# Flanks must be READ-supported, not merely present in the VCF. A block can
# hold sites and no tagged reads at all -- the read-tagging margin drops reads
# that observe too few of its sites -- and a flank with no reads cannot be
# linked by read identity, however many sites it has. Measured on
# chr20:48,176,830-48,229,446, where the nearest left block (48162480, 2 sites)
# carries zero reads, so the first version of this script reported n=0 shared
# reads and called the gap unclosable.
base_reads_per_ps = collections.Counter(ps for _, ps in base_tags.values())
base_read_span = {}
for q, (h, ps) in base_tags.items():
    lo_hi = base_read_span.setdefault(ps, [10 ** 12, 0])
    e = base_sets.get(ps)
    if e:
        lo_hi[0] = min(lo_hi[0], e[0])
        lo_hi[1] = max(lo_hi[1], e[1])
# The flanks are the blocks holding the phased sites nearest the gap, not the
# blocks with the largest extent or the most reads. Block extents interleave and
# a flank can be a two-site island: at chr20:48,176,830 the adjacent block
# (48162480) holds 2 sites and no tagged reads, while the nearest read-supported
# block is 80 kb further out -- picking by read support silently stitched to the
# wrong side of an 80 kb stretch that holds no phased sites at all. A read-less
# flank is still linkable, through its own genotypes rather than its tags.
supported = {ps for ps, n in base_reads_per_ps.items()
             if n >= a.min_flank_reads and ps in base_sets}
base_vcf_sites = vcf_sites(base_vcf)
left_site = max((pos for pos in base_vcf_sites if pos < GL), default=None)
right_site = min((pos for pos in base_vcf_sites if pos > GR), default=None)
left_ps = base_vcf_sites[left_site][3] if left_site else None
right_ps = base_vcf_sites[right_site][3] if right_site else None

verdict = dict(gap=[GL, GR], gap_bp=GR - GL, region=region,
               baseline_blocks=len(base_sets), gap_blocks=len(gap_sets),
               left_flank_ps=left_ps, right_flank_ps=right_ps)

print(f'gap {GL}-{GR}  ({(GR - GL) / 1e3:.1f} kb)   window {region}')
print(f'  baseline: {len(base_sets)} blocks, {len(base_tags)} tagged reads in '
      f'the gauge window, {len(supported)} with >= {a.min_flank_reads} reads')
print(f'  nearest phased site left of the gap: {left_site} (PS={left_ps}, '
      f'{base_reads_per_ps.get(left_ps, 0)} tagged reads, '
      f'{base_sets.get(left_ps, (0, 0, 0))[2]} sites)')
print(f'  nearest phased site right of the gap: {right_site} (PS={right_ps}, '
      f'{base_reads_per_ps.get(right_ps, 0)} tagged reads, '
      f'{base_sets.get(right_ps, (0, 0, 0))[2]} sites)')
if left_ps:
    e = base_sets[left_ps]
    print(f'    left  {e[0]}-{e[1]} ({e[2]} sites)')
if right_ps:
    e = base_sets[right_ps]
    print(f'    right {e[0]}-{e[1]} ({e[2]} sites)')
spanning_base = [ps for ps, e in base_sets.items() if e[0] <= GL and e[1] >= GR]
verdict['baseline_spans_gap'] = bool(spanning_base)
print(f'  baseline already spans the gap: {"YES" if spanning_base else "no"}')

# ---- 2. EVIDENCE inside the interval.
rows = []
for r in csv.DictReader(open(gapdir / 'candidates.tsv'), delimiter='\t'):
    pos = int(r['POS'])
    if not (GL <= pos <= GR) or 'HET' not in r['CATEGORY']:
        continue
    seg, n = (None, 0)
    if truth:
        seg, n = segregation(pos, r['REF'], r['ALT'], r['TYPE'], a.site_flank, truth)
    rows.append(dict(pos=pos, type=r['TYPE'], category=r['CATEGORY'],
                     phase_set=r['PHASE_SET'], segregation=seg, scored=n))
verdict['gap_het_sites'] = len(rows)
verdict['gap_het_phased'] = sum(1 for r in rows if r['phase_set'] not in ('-1', '', '0', '.'))
informative = [r for r in rows if r['segregation'] is not None and r['segregation'] >= 0.90]
verdict['gap_het_informative'] = len(informative)
print(f'\n  BAM evidence in the interval: {len(rows)} het sites, '
      f'{verdict["gap_het_phased"]} phased by the BAM channel'
      + (f', {len(informative)} informative vs truth' if truth else ''))
print('  %-11s %-6s %-18s %-12s %s' % ('pos', 'type', 'category', 'phase_set', 'segregation'))
for r in rows[:14]:
    seg = f'{r["segregation"]:.3f} (n={r["scored"]})' if r['segregation'] is not None else '--'
    print('  %-11d %-6s %-18s %-12s %s' % (r['pos'], r['type'], r['category'],
                                           r['phase_set'], seg))

# ---- 4. STITCH.
# STITCH, in two stages. The gap channel's own blocks are composed with each
# other first: they overlap in read space, so a read carrying both can orient
# one against the other, and that is the block-to-block link the pipeline has no
# mechanism for (select_stitch_orientation is wired only to the adjacent-chunk
# seam and to gap recovery's flank votes). Only then is each composed frame
# linked to the flanks. Composing second would ask each single block to reach
# both flanks, which is exactly the all-or-nothing test recovery already fails.
# Only blocks whose phased sites lie inside the gap (plus its read margin) are
# the gap's own phasing. A block reaching far outside it is a flank re-solve.
gap_in_window = [ps for ps, e in gap_sets.items()
                 if e[1] >= GL - a.gap_margin and e[0] <= GR + a.gap_margin]
parent = {ps: ps for ps in gap_in_window}
parity = {ps: 0 for ps in gap_in_window}


def find(ps):
    path = []
    flip = 0
    while parent[ps] != ps:
        path.append(ps)
        flip ^= parity[ps]
        ps = parent[ps]
    return ps, flip


print('\n  stage 1 -- compose the gap channel\'s own blocks '
      f'(net-margin rule, margin {a.stitch_margin})')
compositions = []
gap_vcf_sites = vcf_sites(gapdir / 'native.vcf')
bridge_reads = load_reads(max(1, GL - a.seam_read_window),
                          GR + a.seam_read_window)
ordered = sorted(gap_in_window, key=lambda ps: gap_sets[ps][0])
for left, right in zip(ordered, ordered[1:]):
    t, shared = votes(gap_tags, left, gap_tags, right)
    flip = decide(t, a.stitch_margin)
    source = 'tags'
    if flip is None:
        # No read is tagged in both blocks, which is the normal state when two
        # blocks split at the same position. Fall back to the alleles.
        br = bridge_vote({pos: v for pos, v in gap_vcf_sites.items() if v[3] == left},
                         {pos: v for pos, v in gap_vcf_sites.items() if v[3] == right},
                         bridge_reads, a.seam_window, a.seam_sites,
                         a.min_site_obs, a.site_flank, a.seam_span)
        if br is not None:
            t, shared = br['votes'], br['voters']
            flip = decide(t, a.stitch_margin)
            source = f'alleles ({br["crossing"]} reads cross the seam)'
    compositions.append(dict(a=left, b=right, votes=t, shared=shared,
                             flip=flip, source=source))
    print(f'    {left} ({gap_sets[left][2]} sites) + {right} '
          f'({gap_sets[right][2]} sites): n={shared} votes={t} [{source}] -> '
          + ('no link' if flip is None else f'compose, flip={int(flip)}'))
    if flip is None:
        continue
    ra, fa = find(left)
    rb, fb = find(right)
    if ra == rb:
        continue
    parent[rb] = ra
    parity[rb] = fa ^ fb ^ int(flip)

frames = collections.defaultdict(dict)
for ps in gap_in_window:
    root, flip = find(ps)
    frames[root][ps] = flip


def frame_tags(members):
    out = {}
    for ps, flip in members.items():
        for q, (h, p_) in gap_tags.items():
            if p_ == ps:
                out[q] = (3 - h if flip else h, 'frame')
    return out


print(f'  -> {len(frames)} composed frame(s) from {len(gap_in_window)} gap block(s)')
print('\n  stage 2 -- composed frame -> flank votes')
links = []
for root, members in sorted(frames.items(), key=lambda kv: -len(kv[1])):
    ftags = frame_tags(members)
    lo = min(gap_sets[ps][0] for ps in members)
    hi = max(gap_sets[ps][1] for ps in members)
    sites = sum(gap_sets[ps][2] for ps in members)
    row = dict(frame=list(members), span=[lo, hi], sites=sites, reads=len(ftags))
    # The frame's genotypes in its own composed gauge: a member block with
    # parity 1 has its phased genotypes flipped, so the frame reads as one block.
    frame_sites = {}
    for pos, v in gap_vcf_sites.items():
        if v[3] not in members:
            continue
        gt = v[2]
        if members[v[3]]:
            gt = '|'.join(reversed(gt.split('|')))
        frame_sites[pos] = (v[0], v[1], gt, 'frame')
    for side, flank in (('left', left_ps), ('right', right_ps)):
        if flank is None:
            row[side] = None
            continue
        t, shared = votes(ftags, 'frame', base_tags, flank)
        flip = decide(t, a.stitch_margin)
        source = 'tags'
        if flip is None:
            flank_sites = {pos: v for pos, v in base_vcf_sites.items()
                           if v[3] == flank}
            if side == 'left':
                br = bridge_vote(flank_sites, frame_sites, bridge_reads,
                                 a.seam_window, a.seam_sites, a.min_site_obs,
                                 a.site_flank, a.seam_span)
            else:
                br = bridge_vote(frame_sites, flank_sites, bridge_reads,
                                 a.seam_window, a.seam_sites, a.min_site_obs,
                                 a.site_flank, a.seam_span)
            if br is not None:
                t, shared = br['votes'], br['voters']
                flip = decide(t, a.stitch_margin)
                # The vote is (flank, frame) on the left and (frame, flank) on
                # the right; the orientation applied to the frame is the same
                # either way because the table is symmetric under transpose.
                source = (f'alleles ({br["crossing"]} cross, '
                          f'{br["n_left_sites"]}+{br["n_right_sites"]} sites)')
        row[side] = dict(votes=t, shared=shared, flip=flip, source=source)
        print(f'    frame {sorted(members)} ({sites} sites, {lo}-{hi}, '
              f'{len(ftags)} reads) vs {side} flank {flank}: n={shared} '
              f'votes={t} [{source}] -> '
              + ('no link' if flip is None else f'link, flip={int(flip)}'))
    row['members'] = members
    links.append(row)

closed = [r for r in links
          if r.get('left') and r.get('right')
          and r['left']['flip'] is not None and r['right']['flip'] is not None]
# A frame that reaches only one flank still extends it. That is not a closed gap,
# but it is phased coverage the baseline does not have, and it is gated the same
# way -- so it is measured rather than discarded.
one_sided = [r for r in links if r not in closed and
             ((r.get('left') and r['left']['flip'] is not None) or
              (r.get('right') and r['right']['flip'] is not None))]
verdict['compositions'] = compositions
verdict['links'] = [{k: v for k, v in r.items() if k != 'members'} for r in links]
verdict['closed'] = bool(closed)

merged = dict(base_tags)
if closed:
    chosen = max(closed, key=lambda r: r['sites'])
    verdict['closing_frame'] = sorted(chosen['members'])
    lflip, rflip = chosen['left']['flip'], chosen['right']['flip']
    # Everything is expressed in the left flank's gauge: the composed frame is
    # flipped onto it if its vote says so, and the right flank follows the
    # composed parity. Reads the baseline already placed in a flank keep their
    # own label -- the frame may only add reads, never relabel phased ones.
    ftags = frame_tags(chosen['members'])
    for q, (h, _) in ftags.items():
        if q in merged and merged[q][1] in (left_ps, right_ps):
            continue
        merged[q] = (3 - h if lflip else h, left_ps)
    right_parity = (lflip != rflip)
    for q, (h, ps) in list(merged.items()):
        if ps == right_ps:
            merged[q] = (3 - h if right_parity else h, left_ps)
    print(f'\n  CLOSED via composed frame {sorted(chosen["members"])} '
          f'(left flip={int(lflip)}, right flip={int(rflip)}, '
          f'right flank parity={int(right_parity)})')
elif one_sided:
    chosen = max(one_sided, key=lambda r: r['reads'])
    side = 'left' if (chosen.get('left') and chosen['left']['flip'] is not None) else 'right'
    flank = left_ps if side == 'left' else right_ps
    flip = chosen[side]['flip']
    verdict['extended_side'] = side
    verdict['extending_frame'] = sorted(chosen['members'])
    ftags = frame_tags(chosen['members'])
    added = 0
    for q, (h, _) in ftags.items():
        if q in merged:
            continue
        merged[q] = (3 - h if flip else h, flank)
        added += 1
    print(f'\n  NOT CLOSED, but the frame {sorted(chosen["members"])} extends the '
          f'{side} flank {flank} (flip={int(flip)}): {added} reads added')
else:
    print('\n  NOT CLOSED: no gap block links either flank')

# ---- 4b. VALIDATE EACH APPLIED LINK against truth, independently of the gate.
if truth:
    applied = []
    if closed:
        applied = [('left', left_ps, chosen['left']['flip']),
                   ('right', right_ps, chosen['right']['flip'])]
    elif one_sided and 'extended_side' in verdict:
        applied = [(verdict['extended_side'],
                    left_ps if verdict['extended_side'] == 'left' else right_ps,
                    chosen[verdict['extended_side']]['flip'])]
    if applied:
        frame_sites_v = {}
        for pos, v in gap_vcf_sites.items():
            if v[3] not in chosen['members']:
                continue
            gt = v[2]
            if chosen['members'][v[3]]:
                gt = '|'.join(reversed(gt.split('|')))
            frame_sites_v[pos] = (v[0], v[1], gt, 'frame')
        f_or, f_n, _ = site_orientation(sorted(frame_sites_v), frame_sites_v,
                                        bridge_reads, truth, a.site_flank)
        print('\n  link validation against truth '
              f'(frame hap1 carries {f_or} over {f_n} sites)')
        verdict['link_validation'] = []
        for side, flank, flip in applied:
            fl_sites = {pos: v for pos, v in base_vcf_sites.items() if v[3] == flank}
            l_or, l_n, _ = site_orientation(sorted(fl_sites), fl_sites,
                                            bridge_reads, truth, a.site_flank)
            if f_or is None or l_or is None:
                ok = None
            else:
                # flip=0 means the frame keeps its polarity relative to the flank,
                # so the two sides must then carry the same truth haplotype.
                ok = ((l_or == f_or) if not flip else (l_or != f_or))
            verdict['link_validation'].append(
                dict(side=side, flank=flank, flip=bool(flip), flank_hap1=l_or,
                     flank_sites=l_n, frame_hap1=f_or, frame_sites=f_n, correct=ok))
            print(f'    {side} flank {flank}: hap1 carries {l_or} over {l_n} sites, '
                  f'flip={int(flip)} -> '
                  + ('CORRECT' if ok else ('WRONG' if ok is False else 'unscorable')))
        bad = [v for v in verdict['link_validation'] if v['correct'] is False]
        verdict['links_correct'] = not bad

# ---- 5. GATE.
if truth:
    before, after = score(base_tags, truth), score(merged, truth)
    gate = collections.Counter()
    for q in set(before) | set(after):
        gate[(before.get(q), after.get(q))] += 1
    verdict['gate'] = dict(
        tagged_before=len(base_tags), tagged_after=len(merged),
        scored_before=len(before), scored_after=len(after),
        concordant_before=sum(before.values()), concordant_after=sum(after.values()),
        conc_to_disc=gate[(True, False)], disc_to_conc=gate[(False, True)],
        newly_concordant=gate[(None, True)], newly_discordant=gate[(None, False)],
        lost_concordant=gate[(True, None)])
    g = verdict['gate']
    print(f'\n  GATE  tagged {g["tagged_before"]} -> {g["tagged_after"]}   '
          f'concordant {g["concordant_before"]} -> {g["concordant_after"]}')
    print(f'        accuracy {100 * g["concordant_before"] / max(g["scored_before"], 1):.2f}%'
          f' -> {100 * g["concordant_after"] / max(g["scored_after"], 1):.2f}%')
    print(f'        concordant->discordant {g["conc_to_disc"]}   '
          f'newly tagged {g["newly_concordant"]} concordant / '
          f'{g["newly_discordant"]} discordant   lost {g["lost_concordant"]}')
    # Pass means the gate is clean AND something was gained: a closed gap, or
    # phased coverage added by an extension. A clean gate on its own is just the
    # baseline restated.
    passed = (g['conc_to_disc'] == 0 and
              g['concordant_after'] > g['concordant_before'] and
              g['newly_discordant'] <= a.max_new_discordant and
              verdict.get('links_correct', True))
    verdict['pass'] = bool(passed)
    print(f'\n  VERDICT: {"PASS" if passed else "FAIL"}'
          f'  (closed={verdict["closed"]}, '
          f'extended={verdict.get("extended_side", "no")}, '
          f'links_correct={verdict.get("links_correct", "n/a")}, '
          f'gate flips={g["conc_to_disc"]}, '
          f'net concordant {g["concordant_after"] - g["concordant_before"]:+d})')

if a.out:
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w') as f:
        json.dump(dict(verdict=verdict, sites=rows), f, indent=2)
    print(f'\nwrote {a.out}')
