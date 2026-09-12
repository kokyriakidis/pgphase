#!/usr/bin/env python3
"""Phase a VCF using HP/PS tags already present on reads.

Every phasing benchmark in this field (HiPhase, LongPhase, WhatsHap) gives each
tool the *same* variant call set and compares the phasings.  pgphase's graph
pipeline does not work that way -- it phases reads directly from the graph
alignment and emits its own calls -- so a naive comparison scores it on a
different variant set than its competitors, which a reviewer will rightly
object to.

This script removes that objection.  It takes reads already haplotagged by any
means (here, pgphase graph phasing) and transfers those haplotypes onto a
supplied VCF: for each heterozygous record, reads covering the site are grouped
by their HP tag and the alleles they carry, and the record is emitted as 0|1 or
1|0 with a PS drawn from the supporting reads.  The result is a phasing of the
*same* call set the competitors were given, so `whatshap compare` scores all
tools on identical variants.

Usage:
    python3 phase_vcf_from_hp.py <tagged.bam> <in.vcf.gz> <out.vcf> [options]

      --min-reads N     minimum haplotype-informative reads at a site [2]
      --min-ratio F     fraction of a haplotype's reads that must agree [0.7]
      --sample NAME     rename the output sample column
      --region REG      restrict to a region (e.g. chr20)

Sites without enough support are emitted unphased, exactly as a phaser would
leave them, so the comparison is not silently inflated by dropping hard sites.
"""

import argparse
import sys
from collections import Counter, defaultdict

import pysam


def read_allele(pileupread, ref_allele, alt_allele):
    """Which allele does this read carry: 0 (ref), 1 (alt), or None (neither)?

    SNVs are decided by the base at the site.  Indels are decided by the length
    change pysam reports at the preceding aligned position, which is how an
    insertion or deletion is anchored in a pileup.
    """
    if pileupread.is_refskip:
        return None
    dlen = len(alt_allele) - len(ref_allele)
    if dlen == 0:
        if pileupread.is_del or pileupread.query_position is None:
            return None
        base = pileupread.alignment.query_sequence[pileupread.query_position]
        if base == ref_allele[0]:
            return 0
        if base == alt_allele[0]:
            return 1
        return None
    # indel: pysam reports the indel length at the anchor base
    if pileupread.indel == dlen:
        return 1
    if pileupread.indel == 0 and not pileupread.is_del:
        return 0
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--min-ratio", type=float, default=0.7)
    ap.add_argument("--sample")
    ap.add_argument("--region")
    ap.add_argument("--min-mapq", type=int, default=0)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    vin = pysam.VariantFile(args.vcf)
    if args.sample:
        vin.header.samples  # touch before rewriting
    hdr = vin.header.copy()
    if "PS" not in hdr.formats:
        hdr.formats.add("PS", 1, "Integer", "Phase set identifier")
    vout = pysam.VariantFile(args.out, "w", header=hdr)

    stats = Counter()
    region = (args.region,) if args.region else ()
    for rec in vin.fetch(*region):
        stats["records"] += 1
        if len(rec.samples) == 0:
            vout.write(rec); continue
        sample = rec.samples[0]
        gt = sample.get("GT")
        # only heterozygous, biallelic records can be phased
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None \
           or gt[0] == gt[1] or len(rec.alts or ()) != 1:
            vout.write(rec); continue
        stats["het"] += 1

        ref, alt = rec.ref, rec.alts[0]
        # counts[hp][allele]
        counts = defaultdict(Counter)
        ps_votes = Counter()
        for col in bam.pileup(rec.chrom, rec.start, rec.start + 1,
                              truncate=True, min_base_quality=0,
                              stepper="samtools", ignore_overlaps=False):
            if col.reference_pos != rec.start:
                continue
            for pr in col.pileups:
                aln = pr.alignment
                if aln.mapping_quality < args.min_mapq:
                    continue
                if not aln.has_tag("HP"):
                    continue
                hp = aln.get_tag("HP")
                if hp not in (1, 2):
                    continue
                a = read_allele(pr, ref, alt)
                if a is None:
                    continue
                counts[hp][a] += 1
                if aln.has_tag("PS"):
                    ps_votes[int(aln.get_tag("PS"))] += 1

        n1, n2 = sum(counts[1].values()), sum(counts[2].values())
        if n1 < args.min_reads or n2 < args.min_reads:
            stats["unphased_support"] += 1
            vout.write(rec); continue

        # Which allele does each haplotype favour, and how cleanly?
        a1, c1 = counts[1].most_common(1)[0]
        a2, c2 = counts[2].most_common(1)[0]
        if c1 / n1 < args.min_ratio or c2 / n2 < args.min_ratio:
            stats["unphased_mixed"] += 1
            vout.write(rec); continue
        if a1 == a2:
            # both haplotypes favour the same allele -- not a het by these reads
            stats["unphased_same"] += 1
            vout.write(rec); continue

        sample["GT"] = (a1, a2)
        sample.phased = True
        if ps_votes:
            sample["PS"] = ps_votes.most_common(1)[0][0]
        stats["phased"] += 1
        vout.write(rec)

    vout.close()
    het = max(stats["het"], 1)
    print(f"  records={stats['records']:,}  het={stats['het']:,}  "
          f"phased={stats['phased']:,} ({stats['phased']/het*100:.1f}% of het)",
          file=sys.stderr)
    print(f"  left unphased: support={stats['unphased_support']:,} "
          f"mixed={stats['unphased_mixed']:,} same-allele={stats['unphased_same']:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
