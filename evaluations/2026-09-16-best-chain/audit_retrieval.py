#!/usr/bin/env python3
"""Do we retrieve each in-gap site with the right allele, depth and assignment?

Three separable questions, measured per site against the alignment directly:

  allele      does the emitted allele correspond to a real mode of the read
              length distribution, and is every mode with real support emitted?
  depth       does the record's DP match the number of reads that actually
              cover the locus, and do its ref/alt counts match the reads
              carrying those alleles?
  assignment  of the reads the record does count, are they parentally pure --
              i.e. is the allele call itself right, as distinct from complete?

Candidate-table conventions, calibrated against the emitted VCF at shared loci:
the table's POS is the VCF POS + 1, an insertion's ALT holds only the inserted
bases (delta = +len(ALT)), and a deletion's REF holds only the deleted bases
with ALT '.' (delta = -len(REF)).
"""
import argparse
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_search as cs

TRACT_FLANK = 25
MODE_MIN_READS = 5


def delta_of(row):
    if row['TYPE'] == 'INS':
        return len(row['ALT'])
    if row['TYPE'] == 'DEL':
        return -len(row['REF'])
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel', type=Path, required=True)
    p.add_argument('--arm', type=Path, required=True)
    p.add_argument('--bam', type=Path, required=True)
    p.add_argument('--truth-map', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    truth = dict(line.split('\t') for line in a.truth_map.read_text().splitlines())

    with a.panel.open() as stream:
        panel = list(csv.DictReader(stream, delimiter='\t'))
    rows = []
    for w in panel:
        gl, gr = int(w['gap_left']), int(w['gap_right'])
        reads = cs.load_reads(str(a.bam), gl - 60000, gr + 60000, 1)
        with (a.arm / f'w{gl}' / 'candidates.tsv').open() as stream:
            cands = [r for r in csv.DictReader(stream, delimiter='\t')
                     if gl <= int(r['POS']) <= gr]
        # Group candidates by locus so multiallelic loci are visible as such.
        by_locus = collections.defaultdict(list)
        for r in cands:
            by_locus[int(r['POS'])].append(r)

        for pos, recs in sorted(by_locus.items()):
            ref_len = max(len(r['REF']) for r in recs)
            is_snp = all(r['TYPE'] == 'SNP' for r in recs)
            # Only indel records drop the anchor base, so only their POS is the
            # VCF POS + 1; a substitution's POS is the VCF POS itself. Reading a
            # substitution one base to the left scores the wrong column and
            # returns chance-level purity, which is how this was caught.
            if is_snp:
                lo, hi = pos - TRACT_FLANK, pos + TRACT_FLANK
            else:
                lo, hi = pos - 1 - TRACT_FLANK, pos - 1 + ref_len + TRACT_FLANK
            covering = [rd for rd in reads if rd.beg <= lo and rd.end >= hi]
            if len(covering) < 10:
                continue
            obs = {}
            if is_snp:
                for rd in covering:
                    b = rd.base_at(pos)
                    if b is not None:
                        obs[rd.name] = b.upper()
                modes = collections.Counter(obs.values())
                emitted = {r['ALT'].upper() for r in recs} | {recs[0]['REF'].upper()}
            else:
                for rd in covering:
                    obs[rd.name] = rd.net_length(lo, hi)
                modes = collections.Counter(obs.values())
                emitted = {delta_of(r) for r in recs} | {0}
            real_modes = {m for m, n in modes.items() if n >= MODE_MIN_READS}
            unrepresented = sorted(real_modes - emitted,
                                   key=lambda m: -modes[m])
            matched = sum(n for m, n in modes.items() if m in emitted)
            dp = sum(int(r['DP']) for r in recs)
            counted = sum(int(r['REF_COUNT']) + int(r['ALT_COUNT']) for r in recs)
            # Assignment purity: of the reads carrying an emitted ALT, how
            # one-sided is their parental composition?
            pure = []
            for r in recs:
                key = r['ALT'].upper() if is_snp else delta_of(r)
                names = [q for q, v in obs.items() if v == key and q in truth]
                if len(names) >= 5:
                    c = collections.Counter(truth[q] for q in names)
                    pure.append(max(c.values()) / sum(c.values()))
            rows.append(dict(
                gap_left=gl, pos=pos, kind='SNP' if is_snp else 'INDEL',
                records=len(recs), covering=len(covering), dp=dp, counted=counted,
                dp_over_cov=round(dp / len(covering), 3),
                modes_with_support=len(real_modes),
                reads_matching_emitted=matched,
                matched_frac=round(matched / len(covering), 3),
                unrepresented_modes=len(unrepresented),
                biggest_unrepresented=(f'{unrepresented[0]}:{modes[unrepresented[0]]}'
                                       if unrepresented else '-'),
                assignment_purity=round(min(pure), 3) if pure else -1.0,
                category=recs[0]['CATEGORY']))

    with a.output.open('w', newline='') as stream:
        wr = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                            lineterminator='\n')
        wr.writeheader()
        wr.writerows(rows)

    print('%-11s %6s %6s %8s %9s %9s %9s %9s' % (
        'gap_left', 'sites', 'indel', 'DP/cov', 'matched', 'unrepr', 'purity', 'multiallelic'))
    for gl in sorted({r['gap_left'] for r in rows}):
        g = [r for r in rows if r['gap_left'] == gl]
        ind = [r for r in g if r['kind'] == 'INDEL']
        pur = [r['assignment_purity'] for r in g if r['assignment_purity'] >= 0]
        print('%-11d %6d %6d %8.3f %8.1f%% %9d %9s %9d' % (
            gl, len(g), len(ind),
            sum(r['dp'] for r in g) / max(sum(r['covering'] for r in g), 1),
            100 * sum(r['reads_matching_emitted'] for r in g) / max(sum(r['covering'] for r in g), 1),
            sum(r['unrepresented_modes'] for r in g),
            '%.3f' % (sum(pur) / len(pur)) if pur else '-',
            sum(1 for r in g if r['records'] > 1)))
    ind = [r for r in rows if r['kind'] == 'INDEL']
    snp = [r for r in rows if r['kind'] == 'SNP']
    print()
    for label, grp in (('substitutions', snp), ('indels', ind)):
        if not grp:
            continue
        pur = [r['assignment_purity'] for r in grp if r['assignment_purity'] >= 0]
        print('%-14s n=%-4d DP/coverage %.3f   reads matching an emitted allele %.1f%%   '
              'sites with an unrepresented mode %d (%.0f%%)   mean assignment purity %s' % (
                  label, len(grp),
                  sum(r['dp'] for r in grp) / max(sum(r['covering'] for r in grp), 1),
                  100 * sum(r['reads_matching_emitted'] for r in grp) / max(sum(r['covering'] for r in grp), 1),
                  sum(1 for r in grp if r['unrepresented_modes']),
                  100 * sum(1 for r in grp if r['unrepresented_modes']) / len(grp),
                  '%.3f' % (sum(pur) / len(pur)) if pur else '-'))
    print(f'wrote {a.output}')


if __name__ == '__main__':
    main()
