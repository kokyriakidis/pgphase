#!/usr/bin/env python3
"""Verify, on stock defaults, that every site is injected and represented right.

Six independent checks per panel window, each a pass/fail with the offending
loci named. Truth is applied only to decide what a site IS, never to select
sites, so a pass means the pipeline got there on its own evidence.

  1 dropped      a candidate the alignment channel has and the hybrid does not
  2 duplicated   two emitted records at one position
  3 missing      a catalog claim that read truth calls heterozygous, absent
                 from the hybrid's candidates
  4 allele set   a locus where no read carries the reference and the record
                 names one allele, so it cannot express the locus
  5 depth        a record whose DP is under half the reads covering its window
  6 verdict      a record classified homozygous that read truth calls
                 heterozygous, or the reverse
  7 attributes   an injected candidate whose own fields are inconsistent --
                 DP against its counts, AF against its counts, or a strand
                 tally that does not sum to the count it accompanies. The
                 strand fields are load-bearing: the ONT strand-bias screen in
                 classify_graph_only_candidates declines to test when
                 forward_alt + reverse_alt is zero.
"""
import argparse, collections, csv, os, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '2026-09-16-best-chain'))
import chain_search as cs

FLANK = 25


def load_candidates(path):
    if not os.path.exists(path):
        return {}
    return {(r['POS'], r['TYPE'], r['REF'], r['ALT']): r
            for r in csv.DictReader(open(path), delimiter='\t')}


def emitted_positions(path):
    counts = collections.Counter()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith('#'):
                continue
            counts[line.split('\t', 2)[1]] += 1
    return counts


def alt_reads_by_parent(reads, truth, pos, vtype, ref, alt):
    """Parental split of the reads carrying this record's ALT.

    A substitution is read at its own position; an indel by net length across
    its footprint, anchored at the VCF position (candidate POS minus one).
    """
    calls = collections.Counter()
    if vtype == 'SNP':
        for read in reads:
            if read.beg > pos or read.end <= pos or read.name not in truth:
                continue
            base = read.base_at(pos)
            if base is not None and base.upper() == alt.upper():
                calls[truth[read.name]] += 1
    else:
        lo, hi = pos - 1 - FLANK, pos - 1 + len(ref) + FLANK
        delta = len(alt) if vtype == 'INS' else -len(ref)
        for read in reads:
            if read.beg > lo or read.end < hi or read.name not in truth:
                continue
            if read.net_length(lo, hi) == delta:
                calls[truth[read.name]] += 1
    total = sum(calls.values())
    return (max(calls.values()) / total if total else 0.0), total


def window_composition(reads, truth, pos, ref_len):
    """Net-length composition over reads fully covering the footprint."""
    lo, hi = pos - 1 - FLANK, pos - 1 + ref_len + FLANK
    agg = collections.defaultdict(collections.Counter)
    for read in reads:
        if read.beg > lo or read.end < hi or read.name not in truth:
            continue
        agg[read.net_length(lo, hi)][truth[read.name]] += 1
    return agg, sum(sum(v.values()) for v in agg.values())


def covering_reads(reads, pos, ref_len):
    lo, hi = pos - 1 - FLANK, pos - 1 + ref_len + FLANK
    return sum(1 for read in reads if read.beg <= lo and read.end >= hi)


def catalog_claims(path, lo, hi):
    out = []
    proc = subprocess.run(['tabix', path, f'CHM13#0#chr20:{lo}-{hi}'],
                          capture_output=True, text=True)
    for line in proc.stdout.splitlines():
        f = line.split('\t')
        for alt in f[4].split(','):
            out.append((int(f[1]), f[3], alt))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel', required=True)
    p.add_argument('--hybrid', required=True, help='arm directory, stock defaults')
    p.add_argument('--bam-arm', required=True, help='collect-bam-variation arm directory')
    p.add_argument('--catalog', required=True)
    p.add_argument('--bam', required=True)
    p.add_argument('--truth-map', required=True)
    p.add_argument('--min-truth-reads', type=int, default=10)
    p.add_argument('--output', required=True)
    a = p.parse_args()

    truth = dict(l.split('\t') for l in open(a.truth_map).read().splitlines())
    rows = []
    findings = collections.defaultdict(list)

    for w in csv.DictReader(open(a.panel), delimiter='\t'):
        gl, gr = int(w['gap_left']), int(w['gap_right'])
        lo, hi = gl - 50000, gr + 50000
        reads = cs.load_reads(a.bam, lo - 2000, hi + 2000, 1)
        hyb = load_candidates(os.path.join(a.hybrid, f'w{gl}', 'candidates.tsv'))
        bam = load_candidates(os.path.join(a.bam_arm, f'w{gl}', 'candidates.tsv'))

        dropped = []
        for key in set(bam) - set(hyb):
            pos, vtype, ref, alt = key
            purity, n = alt_reads_by_parent(reads, truth, int(pos), vtype, ref, alt)
            if n >= a.min_truth_reads and purity >= 0.90:
                dropped.append((pos, vtype, round(purity, 3), n))

        dup = [pos for pos, n in emitted_positions(
            os.path.join(a.hybrid, f'w{gl}', 'native.vcf')).items() if n > 1]

        have = {(int(k[0]), int(k[0]) - 1) for k in hyb}
        have_pos = {p for pair in have for p in pair}
        missing = []
        for pos, ref, alt in catalog_claims(a.catalog, gl, gr):
            if pos in have_pos:
                continue
            vtype = 'SNP' if len(ref) == 1 and len(alt) == 1 else (
                'INS' if len(alt) > len(ref) else 'DEL')
            purity, n = alt_reads_by_parent(reads, truth, pos, vtype, ref, alt)
            if n >= a.min_truth_reads and purity >= 0.90:
                missing.append((pos, vtype, round(purity, 3), n))

        thin, single, verdict = [], [], []
        for key, r in hyb.items():
            pos, vtype, ref, alt = int(key[0]), key[1], key[2], key[3]
            if not (gl - 2000 <= pos <= gr + 2000):
                continue
            ref_len = len(ref) if vtype == 'DEL' else 1
            cov = covering_reads(reads, pos, ref_len)
            dp = int(r['DP'])
            if cov >= 20 and dp < 0.5 * cov:
                thin.append((key[0], vtype, dp, cov))
            if vtype != 'SNP':
                agg, total = window_composition(reads, truth, pos, ref_len)
                if total >= 20:
                    at_ref = sum(agg.get(0, collections.Counter()).values())
                    # Two clean parental modes are required before this counts as
                    # a representation failure. Without them the locus has no
                    # allele pair to express: in a repeat tract the reads scatter
                    # across many net lengths and neither haplotype has a modal
                    # one, so a single-allele record is not the wrong description
                    # of a two-allele locus, it is the only description available.
                    # Measured over the panel, 48 of 51 loci that carry no
                    # reference read are of that kind, and flagging them buries
                    # the 3 that are real.
                    clean_modes = [d for d, c in agg.items()
                                   if d != 0 and sum(c.values()) >= 8 and
                                   max(c.values()) / sum(c.values()) >= 0.90]
                    if at_ref < 0.10 * total and len(clean_modes) >= 2 and \
                       ',' not in r['ALT'] and not r['ALT'].endswith(','):
                        single.append((key[0], vtype, r['CATEGORY'],
                                       sorted(clean_modes)))
            purity, n = alt_reads_by_parent(reads, truth, pos, vtype, ref, alt)
            if n >= a.min_truth_reads:
                is_het = purity >= 0.90
                called_hom = 'HOM' in r['CATEGORY']
                if is_het and called_hom:
                    verdict.append((key[0], vtype, r['CATEGORY'], round(purity, 3), n))

        attrs = []
        for key in set(hyb) - set(bam):
            r = hyb[key]
            dp, rc, ac = int(r['DP']), int(r['REF_COUNT']), int(r['ALT_COUNT'])
            fr, rr = int(r['FORWARD_REF']), int(r['REVERSE_REF'])
            fa, ra = int(r['FORWARD_ALT']), int(r['REVERSE_ALT'])
            why = []
            if dp != rc + ac:
                why.append('DP != ref+alt')
            if abs(float(r['AF']) - (ac / (rc + ac) if rc + ac else 0.0)) > 1e-6:
                why.append('AF != alt/(ref+alt)')
            if fa + ra != ac:
                why.append('alt strand %d+%d != %d' % (fa, ra, ac))
            if fr + rr != rc:
                why.append('ref strand %d+%d != %d' % (fr, rr, rc))
            if why:
                attrs.append((key[0], key[1], '; '.join(why)))

        for name, lst in (('dropped', dropped), ('duplicated', dup), ('missing', missing),
                          ('allele_set', single), ('depth', thin), ('verdict', verdict),
                          ('attributes', attrs)):
            for item in lst:
                findings[name].append((gl, item))
        rows.append(dict(gap_left=gl, candidates=len(hyb), dropped=len(dropped),
                         duplicated=len(dup), missing=len(missing),
                         allele_set=len(single), depth=len(thin), verdict=len(verdict),
                         attributes=len(attrs)))

    with open(a.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    cols = ('candidates', 'dropped', 'duplicated', 'missing', 'allele_set',
            'depth', 'verdict', 'attributes')
    print('%-11s %10s %8s %10s %8s %10s %6s %8s %11s' % (
        'gap_left', 'candidates', 'dropped', 'duplicated', 'missing',
        'alleleSet', 'depth', 'verdict', 'attributes'))
    for r in rows:
        print('%-11d %10d %8d %10d %8d %10d %6d %8d %11d' % (
            r['gap_left'], *[r[c] for c in cols]))
    print('%-11s %10d %8d %10d %8d %10d %6d %8d %11d' % ('TOTAL',
        *[sum(r[c] for r in rows) for c in cols]))
    print()
    for name in ('dropped', 'duplicated', 'missing', 'allele_set', 'depth',
                 'verdict', 'attributes'):
        items = findings[name]
        print('%s: %d' % (name.upper(), len(items)))
        for gl, item in items[:8]:
            print('   w%-10d %s' % (gl, item))
        if len(items) > 8:
            print('   ... and %d more' % (len(items) - 8))
    print('wrote %s' % a.output)


if __name__ == '__main__':
    main()
