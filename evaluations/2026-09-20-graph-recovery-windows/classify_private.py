#!/usr/bin/env python3
"""Count BAM recovery alleles absent from the raw graph catalog by edited sequence."""

import argparse
import csv

import pysam


def catalog_key(pos, ref, alt):
    if len(ref) == len(alt):
        prefix = 0
        while prefix < len(ref) and ref[prefix] == alt[prefix]:
            prefix += 1
        suffix = 0
        while suffix + prefix < len(ref) and ref[-suffix - 1] == alt[-suffix - 1]:
            suffix += 1
        if len(ref) - prefix - suffix == 1:
            return pos + prefix, 8, 1, alt[prefix:prefix + 1]
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    return pos + prefix, 1 if len(alt) > len(ref) else 2, len(ref) - prefix, alt[prefix:]


def equivalent(a, b, fasta, chrom):
    if a[1] == 8 or b[1] == 8:
        return a == b
    if abs(a[0] - b[0]) > 100 or len(a[3]) - a[2] != len(b[3]) - b[2]:
        return False
    beg = max(1, min(a[0], b[0]) - 1)
    end = max(a[0] + a[2], b[0] + b[2]) + 1
    reference = fasta.fetch(chrom, beg - 1, end).upper()

    def edited(key):
        pos, _, ref_len, alt = key
        offset = pos - beg
        return reference[:offset] + alt + reference[offset + ref_len:]

    return edited(a) == edited(b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audit', help='--recovery-audit-out TSV from collect-graph-variation')
    parser.add_argument('--sites', default='test_data/chr20.sites.striped.vcf.gz')
    parser.add_argument('--ref', default='test_data/chm13v2.0.chr20.renamed.fa')
    parser.add_argument('--chrom', default='CHM13#0#chr20')
    args = parser.parse_args()
    sites = pysam.VariantFile(args.sites)
    fasta = pysam.FastaFile(args.ref)
    represented = []
    private = []
    for row in csv.DictReader(open(args.audit), delimiter='\t'):
        if (row['INSIDE_WINDOW'] != '1' or row['KNOWN_RAW'] == '1' or
                row['KNOWN_TRANSLATED'] == '1' or row['ALN_VERIFIED'] != '1' or
                row['CATEGORY'] not in ('3', '4', '6')):
            continue
        kind = int(row['TYPE'])
        pos = int(row['POS']) + (kind != 8)
        allele = (pos, kind, int(row['REF_LEN']),
                  '' if row['ALT'] == '.' else row['ALT'])
        match = False
        for site in sites.fetch(args.chrom, max(0, pos - 101), pos + 101):
            for alt in site.alts or ():
                if alt == '*' or alt.startswith('<'):
                    continue
                if equivalent(allele, catalog_key(site.pos, site.ref, alt),
                              fasta, args.chrom):
                    match = True
                    break
            if match:
                break
        (represented if match else private).append(row['POS'])
    print(f'verified unmatched: {len(represented) + len(private)}')
    print(f'represented in raw catalog: {len(represented)}')
    print(f'absent from raw catalog: {len(private)}')
    print('absent positions:', ','.join(private) or '.')


if __name__ == '__main__':
    main()
