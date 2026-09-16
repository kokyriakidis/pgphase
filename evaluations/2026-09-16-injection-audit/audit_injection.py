#!/usr/bin/env python3
"""Audit one window's BAM-site injection in four stages, each pass/fail.

The reference set is what `collect-bam-variation` itself calls in the interval:
same reads, same reference, no graph channel. Every stage then asks one question
about that set and fails loudly rather than reporting a number to be read:

  1. PRESENT  -- does the hybrid chunk hold every site the BAM channel calls?
                 A site the BAM channel calls and the chunk lacks is a discovery
                 failure, not an admission one.
  2. FIELDS   -- are depth, allele counts, allele fraction, type, homopolymer
                 flag and category identical between the channels? Any
                 divergence is a transfer bug.
  3. ADMITTED -- is the site's category non-zero when phasing runs? Under
                 --recover-gaps every non-graph candidate is zeroed before the
                 solve, so this is what the retry is supposed to undo, inside
                 the failed window only.
  4. USED     -- did the solve give it a usable heterozygous genotype, do its
                 alleles segregate with read truth, and do the reads it phases
                 keep their tags in the emitted BAM?

Stage 4 carries its own control: clean het classes must come out informative. If
they do not, the genotyping is wrong rather than the sites.
"""
import argparse
import collections
import csv
import json
import re
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bam-candidates', type=Path, required=True,
               help="collect-bam-variation candidates.tsv for the interval")
p.add_argument('--bam-vcf', type=Path, required=True)
p.add_argument('--hybrid-candidates', type=Path, required=True)
p.add_argument('--hybrid-vcf', type=Path, required=True)
p.add_argument('--hybrid-phased-bam', type=Path, required=True)
p.add_argument('--alignment', type=Path, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--left', type=int, required=True)
p.add_argument('--right', type=int, required=True)
p.add_argument('--min-scored', type=int, default=20)
p.add_argument('--flank', type=int, default=25,
               help='bp each side of the site a read must span to be genotyped')
p.add_argument('--min-segregation', type=float, default=0.90)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

HET = ('CLEAN_HET_SNP', 'CLEAN_HET_INDEL', 'NOISY_CAND_HET', 'REP_HET_INDEL')
CLEAN = ('CLEAN_HET_SNP', 'CLEAN_HET_INDEL')


def candidates(path):
    """(pos, type, ref, alt) -> row, restricted to the interval."""
    out = {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            pos = int(r['POS'])
            if not (a.left <= pos <= a.right):
                continue
            out[(pos, r['TYPE'], r['REF'], r['ALT'])] = r
    return out


def vcf_genotypes(path):
    """pos -> (ref, alt, gt, ps) for every record, phased or not."""
    out = collections.defaultdict(list)
    with path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or not (a.left <= int(c[1]) <= a.right):
                continue
            keys, vals = c[8].split(':'), c[9].split(':')
            ps = vals[keys.index('PS')].strip() if 'PS' in keys else None
            out[int(c[1])].append((c[3], c[4], vals[0], ps))
    return out


truth = dict(l.split('\t') for l in a.truth_map.read_text().splitlines())

CONSUME_REF = frozenset('MDN=X')


def reads_over(lo, hi):
    out = subprocess.run(['samtools', 'view', '-q', '1', str(a.alignment),
                          f'{a.contig}:{lo}-{hi}'], capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f'samtools view failed: {out.stderr[:300]}')
    rows = []
    for line in out.stdout.splitlines():
        f = line.split('\t')
        rows.append((f[0], int(f[3]), f[5], f[9]))
    return rows


def net_indel(start, cigar, lo, hi):
    """Net inserted-minus-deleted bases over [lo, hi]; the aligner places a
    repeat-context indel arbitrarily within its run, so it cannot be read at
    the anchor position."""
    ref, delta = start, 0
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(n)
        if op == 'I':
            if lo <= ref <= hi:
                delta += n
        elif op in 'DN':
            if ref + n >= lo and ref <= hi:
                delta -= n
            ref += n
        elif op in CONSUME_REF:
            ref += n
    return delta, start, ref


def base_at(start, cigar, seq, pos):
    ref, q = start, 0
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(n)
        if op in 'M=X':
            if ref <= pos < ref + n:
                return seq[q + pos - ref]
            ref += n
            q += n
        elif op in 'DN':
            if ref <= pos < ref + n:
                return None
            ref += n
        elif op in 'IS':
            q += n
    return None


def segregation(row):
    """How well this site's allele partition follows the read-level truth.

    Genotyping is the validated implementation from
    evaluations/2026-09-16-current-deficit/genotype_sites.py, whose clean-het
    control reads 100% informative: substitutions from the base at the site,
    indels from the net insert-minus-delete length over the site plus a flank,
    and the allele chosen by whichever hypothesis the observed length change is
    CLOSER to, with ties left uncalled. An absolute tolerance instead of the
    closer-hypothesis rule scored the clean-het-indel controls at 0.635 and
    0.533 -- which is what that control is for.
    """
    pos, vtype = int(row['POS']), row['TYPE']
    ref, alt = row['REF'], row['ALT']
    span = max(len(ref), 1)
    flank = a.flank
    rows = reads_over(max(1, pos - flank - span), pos + span + flank)
    calls = {}
    if vtype == 'SNP' and len(ref) == 1 and len(alt) == 1:
        for q, start, cigar, seq in rows:
            b = base_at(start, cigar, seq, pos)
            if b is None:
                continue
            if b.upper() == alt.upper():
                calls[q] = 1
            elif b.upper() == ref.upper():
                calls[q] = 0
    else:
        # This TSV writes a deletion as the deleted bases in REF with ALT '.',
        # and an insertion as the inserted sequence in ALT against one anchor
        # base, so len(ALT) - len(REF) is not the length change.
        first_alt = alt.split(',')[0]
        if vtype == 'DEL' or first_alt == '.':
            expect = -len(ref)
        elif vtype == 'INS':
            expect = len(first_alt) - (len(ref) - 1 if len(ref) > 1 else 0)
        else:
            expect = len(first_alt) - len(ref)
        if expect == 0:
            return None, 0
        lo, hi = pos - flank, pos + span + flank
        for q, start, cigar, seq in rows:
            delta, rbeg, rend = net_indel(start, cigar, lo, hi)
            if rbeg > lo or rend < hi:
                continue
            if abs(delta - expect) < abs(delta):
                calls[q] = 1
            elif abs(delta) < abs(delta - expect):
                calls[q] = 0
    pairs = [(v, truth[q]) for q, v in calls.items() if q in truth]
    if len(pairs) < a.min_scored:
        return None, len(pairs)
    c = collections.Counter(pairs)
    same = c[(0, 'MATERNAL')] + c[(1, 'PATERNAL')]
    flip = c[(0, 'PATERNAL')] + c[(1, 'MATERNAL')]
    return max(same, flip) / len(pairs), len(pairs)


def bam_tags(path):
    out = subprocess.run(['samtools', 'view', str(path)], capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f'samtools view failed: {out.stderr[:300]}')
    tags = {}
    for line in out.stdout.splitlines():
        f = line.split('\t')
        hp = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hp = x[5:]
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hp and ps:
            tags[f[0]] = (int(hp), ps)
    return tags


bam_c, hyb_c = candidates(a.bam_candidates), candidates(a.hybrid_candidates)
# Locus index for presence: (pos, type). Exact allele equality is the wrong test,
# because one event can be represented with a different deletion length on each
# side -- split_nested_msa_deletions rewrites exactly that.
hyb_by_locus = {}
for (pos_, type_, ref_, alt_), row_ in hyb_c.items():
    hyb_by_locus.setdefault((pos_, type_), []).append((ref_, alt_, row_))
bam_gt, hyb_gt = vcf_genotypes(a.bam_vcf), vcf_genotypes(a.hybrid_vcf)
tags = bam_tags(a.hybrid_phased_bam)

FIELDS = ('DP', 'REF_COUNT', 'ALT_COUNT', 'AF', 'TYPE', 'CATEGORY')
present_fields = [f for f in FIELDS if bam_c and f in next(iter(bam_c.values()))]

findings = []
rows = []
# One verdict per LOCUS, not per candidate row. The BAM channel emits both
# nested forms of a repeat deletion at the same position -- that is its own
# defect, and the reason the splitter exists -- so iterating rows double-reports
# every such locus. Where a locus has several representations, the one that
# segregates best is the locus's answer: the question is whether the phase
# information survived injection, not which spelling carried it.
seen_loci = set()
for key, ref_row in sorted(bam_c.items(), key=lambda kv: (kv[0][0], -len(kv[0][2]))):
    if ref_row['CATEGORY'] in HET:
        if (key[0], key[1]) in seen_loci:
            continue
        seen_loci.add((key[0], key[1]))
for key, ref_row in sorted(bam_c.items(), key=lambda kv: (kv[0][0], -len(kv[0][2]))):
    pos, vtype, ref, alt = key
    if ref_row['CATEGORY'] in HET and (pos, vtype, 'done') in seen_loci:
        continue
    if ref_row['CATEGORY'] in HET:
        seen_loci.add((pos, vtype, 'done'))
    hyb_row = hyb_c.get(key)
    same_locus = hyb_by_locus.get((pos, vtype), [])
    representation_differs = hyb_row is None and bool(same_locus)
    if hyb_row is None and same_locus:
        hyb_row = same_locus[0][2]
    rec = dict(pos=pos, type=vtype, ref=ref[:24], alt=alt[:24],
               bam_category=ref_row['CATEGORY'], present=hyb_row is not None,
               hybrid_category=hyb_row['CATEGORY'] if hyb_row else None,
               field_diffs={}, genotype=None, ps=None, vcf_pos=None,
               segregation=None, scored=0, used=False,
               representation_differs=False, locus_rows=[], used_pos=None,
               unscorable_reason=None, unsplit_segregation=None)
    rec['representation_differs'] = representation_differs
    # A nested repeat deletion is legitimately re-represented by
    # split_nested_msa_deletions as the common deletion, homozygous because both
    # haplotypes carry it, plus a residual het a few bases away. The phase
    # information then lives in the RESIDUAL, so a locus must be judged over its
    # whole neighbourhood: keying on the BAM channel's allele reports the common
    # part's 1|1 as a lost site. (At these loci the BAM channel is the one in
    # error -- it emits both nested deletions at one position as independent
    # contradictory hets, which is what the splitter exists to prevent.)
    reach = max(len(ref), 1) + a.flank
    neighbours = [(k, v) for k, v in hyb_c.items()
                  if k[1] == vtype and abs(k[0] - pos) <= reach]
    rec['locus_rows'] = [dict(pos=k[0], ref=k[2][:24], alt=k[3][:24],
                              category=v['CATEGORY']) for k, v in neighbours]
    if hyb_row is None:
        if ref_row['CATEGORY'] in HET:
            findings.append(dict(stage='PRESENT', pos=pos,
                                 detail=f"{ref_row['CATEGORY']} {vtype} called by the BAM channel "
                                        "is absent from the hybrid chunk"))
    else:
        for f in present_fields:
            if ref_row.get(f) != hyb_row.get(f):
                rec['field_diffs'][f] = [ref_row.get(f), hyb_row.get(f)]
        if rec['field_diffs'] and ref_row['CATEGORY'] in HET and not representation_differs:
            findings.append(dict(stage='FIELDS', pos=pos,
                                 detail=f"{ref_row['CATEGORY']} differs between channels: "
                                        f"{rec['field_diffs']}"))
        # An indel's candidate row is keyed on the first deleted or inserted
        # base; its VCF record is keyed on the anchor base before it. Looking
        # the genotype up at the TSV position alone reports every indel as
        # absent, which is what made the clean-het control fail.
        gts = []
        for delta in (0, -1, 1):
            gts = hyb_gt.get(pos + delta, [])
            if gts:
                rec['vcf_pos'] = pos + delta
                break
        match = [g for g in gts if g[0] == ref.replace('.', '') or g[1] == alt] or gts
        if match:
            rec['genotype'], rec['ps'] = match[0][2], match[0][3]
        if ref_row['CATEGORY'] in HET:
            # Which row of this locus actually carries the phase information?
            usable_rows = []
            for (kpos, ktype, kref, kalt), krow in neighbours:
                gts = []
                for delta in (0, -1, 1):
                    gts = hyb_gt.get(kpos + delta, [])
                    if gts:
                        break
                for g in gts:
                    if '|' in g[2] and g[2] not in ('0|0', '1|1') and g[3]:
                        usable_rows.append((kpos, krow, g))
                        break
            rec['used'] = bool(usable_rows)
            if usable_rows:
                kpos, krow, g = usable_rows[0]
                rec['used_pos'], rec['genotype'], rec['ps'] = kpos, g[2], g[3]
                # A split residual CANNOT be scored by an independent net-length
                # test: its window necessarily contains the common deletion the
                # split removed, so a read carrying only the common part shows
                # that length here and is called alt as well. At
                # chr20:48,177,781 the reads carry -8 (27) and -20 (21); testing
                # the 12 bp residual calls both alt, 66 vs 8, and reads 0.500 --
                # an artifact of the measurement, not an uninformative site.
                # Declare it unscorable and judge such a locus on the unsplit
                # allele instead.
                if kpos != pos and representation_differs:
                    rec['unscorable_reason'] = 'residual_window_contains_the_common_deletion'
                    seg, scored = segregation(ref_row)
                    rec['unsplit_segregation'] = seg
                    seg = None
                else:
                    seg, scored = segregation(krow)
            else:
                seg, scored = segregation(ref_row)
            rec['segregation'], rec['scored'] = seg, scored
            if not usable_rows:
                findings.append(dict(stage='USED', pos=pos,
                                     detail=f"{ref_row['CATEGORY']}: no row of this locus is emitted "
                                            f"as a phased het (nearest genotype {rec['genotype']}), "
                                            "so the locus links nothing"))
            elif rec['unscorable_reason']:
                findings.append(dict(stage='UNSCORABLE', pos=rec['used_pos'],
                                     detail=f"used as a phased het; {rec['unscorable_reason']}, so "
                                            f"this audit cannot score it (the unsplit allele reads "
                                            f"{rec['unsplit_segregation']:.3f})"
                                            if rec['unsplit_segregation'] is not None else
                                            f"used as a phased het; {rec['unscorable_reason']}"))
            elif seg is not None and seg < a.min_segregation:
                findings.append(dict(stage='USED', pos=rec['used_pos'],
                                     detail=f"used as a phased het but segregates {seg:.3f} against "
                                            f"read truth over {scored} reads"))
    rows.append(rec)

# Stage 4 control: the classes we already trust must come out informative.
for r in rows:
    if r['bam_category'] in CLEAN and r['segregation'] is None and r['present']:
        seg, scored = segregation(bam_c[(r['pos'], r['type'],
                                         [k for k in bam_c if k[0] == r['pos']][0][2],
                                         [k for k in bam_c if k[0] == r['pos']][0][3])])
        r['segregation'], r['scored'] = seg, scored
control = [r for r in rows if r['bam_category'] in CLEAN and r['segregation'] is not None]
control_ok = bool(control) and all(r['segregation'] >= a.min_segregation for r in control)
if control and not control_ok:
    findings.append(dict(stage='CONTROL', pos=0,
                         detail='clean het controls do not come out informative, so the genotyping '
                                'in this audit is wrong rather than the sites'))

het_rows = [r for r in rows if r['bam_category'] in HET]
verdict = dict(
    window=[a.left, a.right],
    bam_channel_sites=len(bam_c), bam_channel_het=len(het_rows),
    present=sum(1 for r in het_rows if r['present']),
    absent=sum(1 for r in het_rows if not r['present']),
    field_mismatches=sum(1 for r in het_rows if r['field_diffs']),
    used_as_het=sum(1 for r in het_rows if r['used']),
    genotype_hom=sum(1 for r in het_rows if r['present'] and r['genotype'] in ('1|1', '0|0')),
    informative=sum(1 for r in het_rows if (r['segregation'] or 0) >= a.min_segregation),
    unscorable=sum(1 for r in het_rows if r['segregation'] is None and r['present']),
    control_sites=len(control), control_ok=control_ok,
    tagged_reads=len(tags),
    findings=findings,
    passes=not findings,
)
a.output.write_text(json.dumps(dict(verdict=verdict, sites=rows), indent=1))

print('window %d-%d   BAM channel: %d sites, %d het' % (
    a.left, a.right, verdict['bam_channel_sites'], verdict['bam_channel_het']))
print('  1 PRESENT   %d of %d het sites in the hybrid chunk (%d absent)' % (
    verdict['present'], verdict['bam_channel_het'], verdict['absent']))
print('  2 FIELDS    %d het sites differ between channels' % verdict['field_mismatches'])
print('  3/4 USED    %d used as a phased het, %d emitted homozygous, %d unscorable' % (
    verdict['used_as_het'], verdict['genotype_hom'], verdict['unscorable']))
print('      control: %d clean het sites, informative: %s' % (
    verdict['control_sites'], 'YES' if control_ok else 'NO'))
for r in rows:
    if r['bam_category'] in CLEAN:
        print('        control %d %-16s seg=%s scored=%d gt=%s' % (
            r['pos'], r['bam_category'],
            ('%.3f' % r['segregation']) if r['segregation'] is not None else 'None',
            r['scored'], r['genotype']))
print('      informative among used: %d' % verdict['informative'])
for f in findings[:14]:
    print('   %-8s %-10s %s' % (f['stage'], f['pos'], f['detail']))
if len(findings) > 14:
    print('   ... %d more findings' % (len(findings) - 14))
print('VERDICT: %s' % ('PASS' if verdict['passes'] else 'FAIL'))
print('wrote %s' % a.output)
