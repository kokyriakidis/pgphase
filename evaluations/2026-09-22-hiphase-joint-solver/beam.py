#!/usr/bin/env python3
"""Bounded best-first prototype for the joint binary MEC gap model."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

MODULE = Path(__file__).with_name("analyze.py")
spec = importlib.util.spec_from_file_location("joint_analyze", MODULE)
joint = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = joint
spec.loader.exec_module(joint)


def solve_beam(variants, observations, sites, left_sites, right_sites,
               flip_right, queue_size, snp_first):
    site_order = sorted(sites, key=lambda site: (variants[site].pos, site))
    rows = {
        name: {site: allele for site, allele in row.items() if site in sites}
        for name, row in observations.items()
    }
    rows = {name: row for name, row in rows.items() if len(row) >= 2}
    read_order = sorted(rows)
    read_index = {name: i for i, name in enumerate(read_order)}

    # A multiplier larger than the total secondary cost makes scalar cost
    # comparison equivalent to minimizing (SNP cost, other-site cost). Indels
    # can bridge a SNP-disconnected gap but cannot overturn a supported SNP link.
    secondary_cost_bound = sum(
        max(1, variants[site].weight)
        for row in rows.values() for site in row
        if variants[site].variant_type != "X"
    )
    snp_multiplier = secondary_cost_bound + 1 if snp_first else 1

    by_site = {}
    for site in site_order:
        entries = []
        weight = max(1, variants[site].weight)
        if variants[site].variant_type == "X":
            weight *= snp_multiplier
        for name, row in rows.items():
            if site in row:
                entries.append((read_index[name], row[site], weight))
        by_site[site] = entries

    mismatches = np.zeros((1, len(read_order)), dtype=np.int32)
    totals = np.zeros(len(read_order), dtype=np.int32)
    bits = [0]
    for depth, site in enumerate(site_order):
        entries = by_site[site]
        indices = np.fromiter((e[0] for e in entries), dtype=np.int64)
        alleles = np.fromiter((e[1] for e in entries), dtype=np.int8)
        weights = np.fromiter((e[2] for e in entries), dtype=np.int32)
        if indices.size:
            totals[indices] += weights

        if site in left_sites:
            choices = (variants[site].hap1,)
        elif site in right_sites:
            choices = (variants[site].hap1 ^ int(flip_right),)
        else:
            choices = (0, 1)

        children = []
        child_bits = []
        for orientation in choices:
            child = mismatches.copy()
            if indices.size:
                disagreement = alleles != orientation
                if np.any(disagreement):
                    child[:, indices[disagreement]] += weights[disagreement]
            children.append(child)
            child_bits.extend(bit | (orientation << depth) for bit in bits)
        mismatches = np.concatenate(children, axis=0)
        bits = child_bits
        costs = np.minimum(mismatches, totals - mismatches).sum(axis=1)
        if costs.size > queue_size:
            keep = np.argpartition(costs, queue_size - 1)[:queue_size]
            order = keep[np.argsort(costs[keep], kind="stable")]
            mismatches = mismatches[order]
            bits = [bits[int(i)] for i in order]

    costs = np.minimum(mismatches, totals - mismatches).sum(axis=1)
    best = int(np.argmin(costs))
    read_costs = {
        name: int(min(mismatches[best, i], totals[i] - mismatches[best, i]))
        for i, name in enumerate(read_order)
    }
    return int(costs[best]), read_costs, len(site_order), len(read_order)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", type=Path)
    parser.add_argument("gap_left", type=int)
    parser.add_argument("gap_right", type=int)
    parser.add_argument("--queue-size", type=int, default=1000)
    parser.add_argument("--fold-mod", type=int, default=0)
    parser.add_argument("--fold-rem", type=int, default=0)
    parser.add_argument("--snp-first", action="store_true")
    parser.add_argument("--trusted-sites", action="store_true")
    parser.add_argument("--af-margin", type=float, default=0.15)
    parser.add_argument("--full-flanks", action="store_true")
    args = parser.parse_args()

    variants, reads, observations = joint.load_matrix(args.matrix)
    if args.fold_mod:
        def stable_hash(name):
            value = 14695981039346656037
            for byte in name.encode():
                value ^= byte
                value = (value * 1099511628211) & ((1 << 64) - 1)
            return value
        observations = {
            name: row for name, row in observations.items()
            if stable_hash(name) % args.fold_mod == args.fold_rem
        }
    oriented = [
        variant for variant in variants.values()
        if variant.phase_set > 0 and variant.hap1 in (0, 1) and
        variant.hap2 in (0, 1) and variant.hap1 != variant.hap2
    ]
    # Candidate keys can differ by a few bases from VCF coordinates after
    # indel normalization, so first identify the adjacent block labels using
    # the nearest oriented rows.
    left_boundary = min(
        oriented, key=lambda variant: (abs(variant.pos - args.gap_left), -variant.pos)
    )
    right_boundary = min(
        oriented, key=lambda variant: (abs(variant.pos - args.gap_right), variant.pos)
    )
    left_ps, right_ps = left_boundary.phase_set, right_boundary.phase_set

    if args.trusted_sites:
        def centered(variant):
            return abs(variant.allele_fraction - 0.5) <= args.af_margin

        # SNPs are the primary anchors. A centered indel substitutes only when
        # a flank has no usable SNP; additional indels are admitted only if the
        # SNP-only observation graph cannot connect the two anchors.
        left_pool = [v for v in oriented if v.phase_set == left_ps and
                     v.variant_type == "X" and centered(v)]
        right_pool = [v for v in oriented if v.phase_set == right_ps and
                      v.variant_type == "X" and centered(v)]
        if not left_pool:
            left_pool = [v for v in oriented if v.phase_set == left_ps and centered(v)]
        if not right_pool:
            right_pool = [v for v in oriented if v.phase_set == right_ps and centered(v)]
        if not left_pool or not right_pool:
            print("NO_TRUSTED_FLANK")
            return
        left_anchor = min(
            left_pool, key=lambda v: (abs(v.pos - args.gap_left), -v.pos)
        )
        right_anchor = min(
            right_pool, key=lambda v: (abs(v.pos - args.gap_right), v.pos)
        )
        if args.full_flanks:
            seam_beg = min(v.pos for v in left_pool + right_pool)
            seam_end = max(v.pos for v in left_pool + right_pool)
        else:
            seam_beg = min(left_anchor.pos, right_anchor.pos)
            seam_end = max(left_anchor.pos, right_anchor.pos)

        def unique_subset(include_indels):
            selected = {}
            for site, variant in variants.items():
                if not (seam_beg <= variant.pos <= seam_end) or not centered(variant):
                    continue
                if not include_indels and variant.variant_type != "X":
                    continue
                key = (variant.pos, variant.variant_type, variant.ref_len, variant.alt)
                selected.setdefault(key, site)
            return {site: variants[site] for site in selected.values()}

        seam_variants = unique_subset(False)
        retained, left_sites, right_sites = joint.retained_component(
            seam_variants, observations, left_ps, right_ps, True
        )
        if not retained or not left_sites or not right_sites:
            seam_variants = unique_subset(True)
            retained, left_sites, right_sites = joint.retained_component(
                seam_variants, observations, left_ps, right_ps, True
            )
    else:
        left_anchor, right_anchor = left_boundary, right_boundary
        seam_beg = min(left_anchor.pos, right_anchor.pos)
        seam_end = max(left_anchor.pos, right_anchor.pos)
        seam_variants = {
            site: variant for site, variant in variants.items()
            if seam_beg <= variant.pos <= seam_end
        }
        retained, left_sites, right_sites = joint.retained_component(
            seam_variants, observations, left_ps, right_ps, True
        )
    if not retained or not left_sites or not right_sites:
        print("DISCONNECTED")
        return

    results = []
    for flip in (False, True):
        result = solve_beam(
            variants, observations, retained, left_sites, right_sites,
            flip, args.queue_size, args.snp_first,
        )
        results.append(result)
        print(f"flip={int(flip)} cost={result[0]} sites={result[2]} reads={result[3]}")
    if results[0][0] == results[1][0]:
        print("TIE")
        return
    predicted_flip = results[1][0] < results[0][0]
    shared = results[0][1].keys() & results[1][1].keys()
    same_votes = sum(results[0][1][name] < results[1][1][name] for name in shared)
    flip_votes = sum(results[1][1][name] < results[0][1][name] for name in shared)
    print(
        f"WIN flip={int(predicted_flip)} delta={abs(results[0][0]-results[1][0])} "
        f"votes={same_votes}/{flip_votes} vote_agree={((flip_votes > same_votes) == predicted_flip) and same_votes != flip_votes}"
    )


if __name__ == "__main__":
    main()
