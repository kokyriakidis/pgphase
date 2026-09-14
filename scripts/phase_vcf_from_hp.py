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
      --sample NAME     sample to phase (input must be single-sample)
      --region REG      restrict to a region (e.g. chr20)

Sites without enough support are emitted unphased, exactly as a phaser would
leave them, so the comparison is not silently inflated by dropping hard sites.

One subtlety drives the implementation.  An HP tag is only meaningful *within*
a phase set: HP=1 in one phase set and HP=1 in the next are unrelated, since
each block's haplotype labelling is arbitrary.  Pooling reads by HP alone
therefore scrambles every site near a block boundary, where reads from two
phase sets cover the same position with opposite polarity -- which shows up as
sites that look "mixed" and as phase sets that appear to interleave.  Reads are
grouped by (PS, HP) and the site is decided by the phase set contributing the
most haplotype-informative reads.
"""

import argparse
import bisect
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam
from allele_observations import observe_allele


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


def resolve_bam_contig(vcf_contig, bam_references):
    """Resolve chr20 against pangenome-style BAM names such as CHM13#0#chr20."""
    if vcf_contig in bam_references:
        return vcf_contig
    matches = [name for name in bam_references if name.split("#")[-1] == vcf_contig]
    return matches[0] if len(matches) == 1 else None


def variant_key(rec):
    return rec.chrom, rec.start, rec.ref, ",".join(rec.alts)


def vcf_fetch_args(region):
    if not region:
        return ()
    if ":" not in region:
        return (region,)
    contig, interval = region.rsplit(":", 1)
    if "-" not in interval:
        return (contig, int(interval) - 1, int(interval))
    start, end = interval.replace(",", "").split("-", 1)
    return contig, int(start) - 1, int(end)


def alignment_allele(aln, ref_pos, ref_allele, alt_allele):
    """Return the pileup-equivalent allele at one zero-based reference position."""
    query_pos = 0
    reference_pos = aln.reference_start
    cigartuples = aln.cigartuples or ()
    for op_i, (op, length) in enumerate(cigartuples):
        if op in (0, 7, 8):  # M, =, X
            segment_end = reference_pos + length
            if reference_pos <= ref_pos < segment_end:
                offset = ref_pos - reference_pos
                site_query_pos = query_pos + offset
                dlen = len(alt_allele) - len(ref_allele)
                if dlen == 0:
                    seq = aln.query_sequence
                    if seq is None or site_query_pos >= len(seq):
                        return None
                    base = seq[site_query_pos].upper()
                    if base == ref_allele[0].upper():
                        return 0
                    if base == alt_allele[0].upper():
                        return 1
                    return None

                observed_indel = 0
                if ref_pos == segment_end - 1 and op_i + 1 < len(cigartuples):
                    next_op, next_length = cigartuples[op_i + 1]
                    if next_op == 1:  # insertion after the anchor
                        observed_indel = next_length
                    elif next_op == 2:  # deletion after the anchor
                        observed_indel = -next_length
                if observed_indel == dlen:
                    return 1
                if observed_indel == 0:
                    return 0
                return None
            query_pos += length
            reference_pos = segment_end
        elif op in (1, 4):  # insertion, soft clip
            query_pos += length
        elif op in (2, 3):  # deletion, reference skip
            if reference_pos <= ref_pos < reference_pos + length:
                return None
            reference_pos += length
        elif op == 5:  # hard clip
            continue
    return None


def collect_support_by_alignment_scan(bam, vin, region, min_mapq):
    targets = defaultdict(lambda: defaultdict(list))
    # Cache votes use slots in the sorted diploid GT, preserving the original
    # allele indices on output (including 1/2 genotypes).
    genotype_indices = {}
    fetch_region = vcf_fetch_args(region)
    for rec in vin.fetch(*fetch_region):
        if len(rec.samples) == 0:
            continue
        gt = rec.samples[0].get("GT")
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None \
           or gt[0] == gt[1]:
            continue
        key = variant_key(rec)
        genotype_indices[key] = tuple(sorted(gt))
        targets[rec.chrom][rec.start].append(key)

    support = defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))
    for vcf_contig, by_pos in targets.items():
        bam_contig = resolve_bam_contig(vcf_contig, bam.references)
        if bam_contig is None:
            continue
        positions = sorted(by_pos)
        for aln in bam.fetch(bam_contig):
            if aln.is_unmapped or aln.is_secondary or aln.is_qcfail or aln.is_duplicate:
                continue
            if aln.mapping_quality < min_mapq or not aln.has_tag("HP"):
                continue
            hp = aln.get_tag("HP")
            if hp not in (1, 2) or aln.reference_end is None:
                continue
            ps = int(aln.get_tag("PS")) if aln.has_tag("PS") else 0
            left = bisect.bisect_left(positions, aln.reference_start)
            right = bisect.bisect_left(positions, aln.reference_end)
            for pos in positions[left:right]:
                for key in by_pos[pos]:
                    if "," in key[3]:
                        observed = observe_allele(aln, pos, key[2], (key[2], *key[3].split(",")))
                        allele = genotype_indices[key].index(observed) if observed in genotype_indices[key] else None
                    else:
                        allele = alignment_allele(aln, pos, key[2], key[3])
                    if allele is not None:
                        support[key][ps][hp][allele] += 1
    return support


def write_support_cache(path, support):
    # Legacy REF_COUNT/ALT_COUNT columns are GT-slot 0/1 counts for multiallelic
    # records. Rebuild older caches to include their previously omitted records.
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(("CHROM", "POS0", "REF", "ALT", "PS", "HP", "REF_COUNT", "ALT_COUNT"))
        for (chrom, pos, ref, alt), by_ps in support.items():
            for ps, by_hp in by_ps.items():
                for hp, counts in by_hp.items():
                    writer.writerow((chrom, pos, ref, alt, ps, hp,
                                     counts[0], counts[1]))


def read_support_cache(path):
    support = defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))
    with Path(path).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["CHROM"], int(row["POS0"]), row["REF"], row["ALT"])
            counts = support[key][int(row["PS"])][int(row["HP"])]
            counts[0] = int(row["REF_COUNT"])
            counts[1] = int(row["ALT_COUNT"])
    return support


def phase_call(counts, min_reads, min_ratio):
    n1, n2 = sum(counts[1].values()), sum(counts[2].values())
    if n1 < min_reads or n2 < min_reads:
        return None
    a1, c1 = counts[1].most_common(1)[0]
    a2, c2 = counts[2].most_common(1)[0]
    if c1 / n1 < min_ratio or c2 / n2 < min_ratio or a1 == a2:
        return None
    return a1, a2, n1 + n2


def phase_orientation_call(counts, min_reads, min_ratio, allow_single_hap):
    call = phase_call(counts, min_reads, min_ratio)
    if call is not None or not allow_single_hap:
        return call
    single_hap_calls = []
    for hp in (1, 2):
        total = sum(counts[hp].values())
        if total < min_reads:
            continue
        allele, allele_count = counts[hp].most_common(1)[0]
        if allele_count / total < min_ratio:
            continue
        a1 = allele if hp == 1 else 1 - allele
        single_hap_calls.append((a1, 1 - a1, total))
    return single_hap_calls[0] if len(single_hap_calls) == 1 else None


class ParityUnionFind:
    def __init__(self):
        self.parent = {}
        self.parity = {}

    def add(self, item):
        if item not in self.parent:
            self.parent[item] = item
            self.parity[item] = 0

    def find(self, item):
        self.add(item)
        if self.parent[item] != item:
            parent, parent_parity = self.find(self.parent[item])
            self.parity[item] ^= parent_parity
            self.parent[item] = parent
        return self.parent[item], self.parity[item]

    def union(self, left, right, relation):
        left_root, left_parity = self.find(left)
        right_root, right_parity = self.find(right)
        if left_root == right_root:
            return (left_parity ^ right_parity) == relation
        self.parent[right_root] = left_root
        self.parity[right_root] = left_parity ^ right_parity ^ relation
        return True


def build_phase_set_merges(support, min_reads, min_ratio,
                           min_sites, min_margin, allow_single_hap,
                           require_full_side, min_read_support, edge_order,
                           long_edge_distance, long_edge_min_read_support,
                           excluded_edges):
    edge_votes = defaultdict(Counter)
    edge_read_support = defaultdict(Counter)
    edge_positions = defaultdict(lambda: defaultdict(list))
    for key, by_ps in support.items():
        calls = []
        for ps, counts in by_ps.items():
            if ps == 0:
                continue
            full_call = phase_call(counts, min_reads, min_ratio)
            call = full_call or phase_orientation_call(
                counts, min_reads, min_ratio, allow_single_hap)
            if call is not None:
                calls.append((ps, call[0], full_call is not None, call[2]))
        for i, (left_ps, left_a1, left_full, left_support) in enumerate(calls):
            for right_ps, right_a1, right_full, right_support in calls[i + 1:]:
                if require_full_side and not (left_full or right_full):
                    continue
                pair_left, pair_right = left_ps, right_ps
                pair_left_a1, pair_right_a1 = left_a1, right_a1
                if pair_left > pair_right:
                    pair_left, pair_right = pair_right, pair_left
                    pair_left_a1, pair_right_a1 = pair_right_a1, pair_left_a1
                pair = (pair_left, pair_right)
                relation = pair_left_a1 ^ pair_right_a1
                edge_votes[pair][relation] += 1
                edge_read_support[pair][relation] += min(left_support, right_support)
                edge_positions[pair][relation].append(key[1] + 1)

    accepted = []
    for pair, votes in edge_votes.items():
        winner, winner_sites = votes.most_common(1)[0]
        loser_sites = votes[1 - winner]
        margin = winner_sites - loser_sites
        winner_read_support = edge_read_support[pair][winner]
        long_edge_ok = long_edge_distance <= 0 \
            or abs(pair[1] - pair[0]) <= long_edge_distance \
            or winner_read_support >= long_edge_min_read_support
        if winner_sites >= min_sites and margin >= min_margin \
                and winner_read_support >= min_read_support \
                and long_edge_ok \
                and pair not in excluded_edges:
            accepted.append((winner_sites, margin, winner_read_support, pair, winner))
    if edge_order == "reads":
        accepted.sort(key=lambda edge: (edge[2], edge[0], edge[1]), reverse=True)
    else:
        accepted.sort(reverse=True)

    accepted_pairs = {edge[3] for edge in accepted}
    merged = ParityUnionFind()
    conflicts = 0
    for _, _, _, (left_ps, right_ps), relation in accepted:
        if not merged.union(left_ps, right_ps, relation):
            conflicts += 1
    report = []
    for pair, votes in edge_votes.items():
        winner, winner_sites = votes.most_common(1)[0]
        report.append((
            pair[0], pair[1], winner, winner_sites, votes[1 - winner],
            edge_read_support[pair][winner],
            ",".join(str(pos) for pos in edge_positions[pair][winner]),
            pair in accepted_pairs,
        ))
    return merged, len(edge_votes), len(accepted), conflicts, report


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--min-ratio", type=float, default=0.7)
    ap.add_argument("--sample", help="Sample to phase (single-sample VCF only).")
    ap.add_argument("--region")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument(
        "--support-cache",
        help=("Reusable TSV of per-site (PS,HP,allele) counts. The first run "
              "builds it with one BAM scan; later threshold sweeps read it."),
    )
    ap.add_argument(
        "--rebuild-support-cache", action="store_true",
        help="Recompute --support-cache even when the TSV already exists.",
    )
    ap.add_argument(
        "--merge-phase-sets", action="store_true",
        help="Merge PS blocks whose overlapping site calls agree in orientation.",
    )
    ap.add_argument(
        "--merge-min-sites", type=int, default=3,
        help="Minimum agreeing overlap sites for a PS merge [3].",
    )
    ap.add_argument(
        "--merge-min-margin", type=int, default=2,
        help="Minimum agree-minus-conflict site margin for a PS merge [2].",
    )
    ap.add_argument(
        "--merge-call-min-reads", type=int,
        help="Per-haplotype read minimum for overlap-edge calls [output setting].",
    )
    ap.add_argument(
        "--merge-call-min-ratio", type=float,
        help="Allele-purity minimum for overlap-edge calls [output setting].",
    )
    ap.add_argument(
        "--merge-min-read-support", type=int, default=0,
        help="Minimum summed winning read support for a PS overlap edge [0].",
    )
    ap.add_argument(
        "--merge-long-edge-distance", type=int, default=0,
        help=("PS-label distance beyond which stronger read support is required; "
              "0 disables the distance gate [0]."),
    )
    ap.add_argument(
        "--merge-long-edge-min-read-support", type=int, default=0,
        help=("Minimum winning read support beyond --merge-long-edge-distance "
              "[0]."),
    )
    ap.add_argument(
        "--merge-edge-order", choices=("sites", "reads"), default="sites",
        help="Priority for resolving parity cycles in the PS graph [sites].",
    )
    ap.add_argument(
        "--merge-exclude-edge", action="append", default=[], metavar="PS1,PS2",
        help="Diagnostic: exclude one PS edge from merging. May be repeated.",
    )
    ap.add_argument(
        "--merge-edge-report",
        help="Write all observed PS overlap edges and support to a TSV.",
    )
    ap.add_argument(
        "--merge-single-hap", action="store_true",
        help=("Allow a clean HP1-only or HP2-only site call to orient PS overlap "
              "edges; use with a multi-site merge threshold."),
    )
    ap.add_argument(
        "--merge-require-full-side", action="store_true",
        help=("For single-haplotype PS edges, require the other PS call at the "
              "bridge site to have both HP1 and HP2 support."),
    )
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    vin = pysam.VariantFile(args.vcf)
    if len(vin.header.samples) != 1:
        ap.error("input VCF must contain exactly one sample")
    input_sample = next(iter(vin.header.samples))
    if args.sample and args.sample != input_sample:
        ap.error(f"sample {args.sample!r} is not the input sample {input_sample!r}")
    hdr = vin.header.copy()
    if "PS" not in hdr.formats:
        hdr.formats.add("PS", 1, "Integer", "Phase set identifier")
    vout = pysam.VariantFile(args.out, "w", header=hdr)

    stats = Counter()
    bam_contigs = {
        contig: resolve_bam_contig(contig, bam.references)
        for contig in vin.header.contigs
    }
    support_cache = None
    if args.support_cache:
        cache_path = Path(args.support_cache)
        if cache_path.exists() and not args.rebuild_support_cache:
            support_cache = read_support_cache(cache_path)
        else:
            support_cache = collect_support_by_alignment_scan(
                bam, vin, args.region, args.min_mapq)
            write_support_cache(cache_path, support_cache)
    phase_set_merges = None
    if args.merge_phase_sets:
        if support_cache is None:
            ap.error("--merge-phase-sets requires --support-cache")
        merge_call_min_reads = args.merge_call_min_reads \
            if args.merge_call_min_reads is not None else args.min_reads
        merge_call_min_ratio = args.merge_call_min_ratio \
            if args.merge_call_min_ratio is not None else args.min_ratio
        excluded_edges = set()
        for spec in args.merge_exclude_edge:
            try:
                left, right = (int(value) for value in spec.split(",", 1))
            except ValueError:
                ap.error(f"invalid --merge-exclude-edge value: {spec}")
            excluded_edges.add(tuple(sorted((left, right))))
        phase_set_merges, n_edges, n_accepted, n_conflicts, edge_report = \
            build_phase_set_merges(
            support_cache, merge_call_min_reads, merge_call_min_ratio,
            args.merge_min_sites, args.merge_min_margin, args.merge_single_hap,
            args.merge_require_full_side, args.merge_min_read_support,
            args.merge_edge_order, args.merge_long_edge_distance,
            args.merge_long_edge_min_read_support, excluded_edges)
        if args.merge_edge_report:
            with Path(args.merge_edge_report).open("w", newline="") as report_fh:
                writer = csv.writer(report_fh, delimiter="\t")
                writer.writerow(("PS1", "PS2", "RELATION", "WIN_SITES",
                                 "LOSE_SITES", "READ_SUPPORT", "POSITIONS",
                                 "ACCEPTED"))
                writer.writerows(edge_report)
        print(f"  PS overlap graph: edges={n_edges:,} accepted={n_accepted:,} "
              f"cycle-conflicts={n_conflicts:,}", file=sys.stderr)
    region = vcf_fetch_args(args.region)
    for rec in vin.fetch(*region):
        stats["records"] += 1
        if len(rec.samples) == 0:
            vout.write(rec); continue
        sample = rec.samples[0]
        gt = sample.get("GT")
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None \
           or gt[0] == gt[1]:
            vout.write(rec); continue
        stats["het"] += 1

        # The input caller may already have phased this genotype. Every output
        # phase must come from the supplied HP/PS tags, including at sites that
        # later fail support or representation checks.
        sample["GT"] = tuple(gt)
        sample.phased = False
        if "PS" in hdr.formats:
            sample["PS"] = None
        genotype_indices = tuple(sorted(gt))
        if any(not set(rec.alleles[i].upper()) <= set("ACGT") for i in genotype_indices):
            stats["unphased_unsupported_allele"] += 1
            vout.write(rec); continue
        stats["eligible_het"] += 1

        ref, alt = rec.ref, rec.alts[0]
        bam_contig = bam_contigs.get(rec.chrom)
        if bam_contig is None:
            stats["unphased_contig"] += 1
            vout.write(rec); continue
        # by_ps[phase_set][hp][allele] -- HP is only comparable within a PS
        if support_cache is not None:
            by_ps = support_cache.get(variant_key(rec), {})
        else:
            by_ps = defaultdict(lambda: defaultdict(Counter))
            for col in bam.pileup(bam_contig, rec.start, rec.start + 1,
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
                    if len(rec.alts) > 1:
                        observed = observe_allele(aln, rec.start, rec.ref, rec.alleles)
                        a = genotype_indices.index(observed) if observed in genotype_indices else None
                    else:
                        a = read_allele(pr, ref, alt)
                    if a is None:
                        continue
                    ps = int(aln.get_tag("PS")) if aln.has_tag("PS") else 0
                    by_ps[ps][hp][a] += 1

        if not by_ps:
            stats["unphased_support"] += 1
            vout.write(rec); continue
        # Decide the site within one phase set: the one with the most
        # haplotype-informative reads here.  Reads from other phase sets carry
        # an unrelated HP polarity and must not vote.
        phase_set, counts = max(
            by_ps.items(), key=lambda kv: sum(sum(c.values()) for c in kv[1].values()))
        if len(by_ps) > 1:
            stats["boundary_sites"] += 1

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

        if phase_set_merges is not None and phase_set:
            phase_set, parity = phase_set_merges.find(phase_set)
            if parity:
                a1, a2 = a2, a1

        sample["GT"] = (genotype_indices[a1], genotype_indices[a2])
        sample.phased = True
        if phase_set:
            sample["PS"] = phase_set
        stats["phased"] += 1
        vout.write(rec)

    vout.close()
    het = max(stats["het"], 1)
    print(f"  records={stats['records']:,}  het={stats['het']:,}  "
          f"eligible-het={stats['eligible_het']:,}  phased={stats['phased']:,} "
          f"({stats['phased']/het*100:.1f}% of het)",
          file=sys.stderr)
    print(f"  left unphased: support={stats['unphased_support']:,} "
          f"mixed={stats['unphased_mixed']:,} same-allele={stats['unphased_same']:,} "
          f"unsupported-allele={stats['unphased_unsupported_allele']:,}",
          file=sys.stderr)
    if stats["unphased_contig"]:
        print(f"  left unphased: unresolved-contig={stats['unphased_contig']:,}",
              file=sys.stderr)
    print(f"  sites covered by >1 phase set (decided within the dominant one): "
          f"{stats['boundary_sites']:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
