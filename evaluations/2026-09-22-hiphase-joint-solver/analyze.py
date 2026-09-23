#!/usr/bin/env python3
"""Solve one dumped recovery seam with the HiPhase binary MEC objective.

This is an evaluation prototype. It uses pgphase's existing graph/BAM allele
observations directly, performs no sequence realignment, and fixes the two graph
flanks in each possible relative orientation. Internal site and read labels are
optimized jointly by an exact MILP. A tied optimum abstains.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.stats import binomtest
from scipy.sparse import coo_matrix

GERMLINE_MASK = 0x004 | 0x008 | 0x080 | 0x100 | 0x200 | 0x1000


@dataclass
class Variant:
    index: int
    pos: int
    variant_type: str
    ref_len: int
    alt: str
    category: int
    weight: int
    phase_set: int
    hap1: int
    hap2: int
    bam_injected: bool
    ref_cov: int = 0
    alt_cov: int = 0
    allele_fraction: float = 0.5


@dataclass
class Read:
    name: str
    hap: int
    phase_set: int


def load_matrix(path: Path):
    variants: dict[int, Variant] = {}
    reads: dict[str, Read] = {}
    observations: dict[str, dict[int, int]] = defaultdict(dict)
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "VAR":
                variant = Variant(
                    index=int(fields[1]),
                    pos=int(fields[2]),
                    variant_type=fields[3],
                    ref_len=int(fields[6]),
                    alt=fields[7],
                    category=int(fields[4]),
                    weight=int(fields[5]),
                    phase_set=int(fields[8]),
                    hap1=int(fields[9]),
                    hap2=int(fields[10]),
                    bam_injected=bool(int(fields[11])),
                    ref_cov=int(fields[14]) if len(fields) > 16 else 0,
                    alt_cov=int(fields[15]) if len(fields) > 16 else 0,
                    allele_fraction=float(fields[16]) if len(fields) > 16 else 0.5,
                )
                variants[variant.index] = variant
            elif fields[0] == "READ":
                reads[fields[1]] = Read(
                    name=fields[1], hap=int(fields[6]), phase_set=int(fields[7])
                )
            elif fields[0] == "OBS":
                allele = int(fields[3])
                if allele in (0, 1):
                    observations[fields[1]][int(fields[2])] = allele
    return variants, reads, observations


def choose_flanks(variants: dict[int, Variant], gap_left: int, gap_right: int):
    anchors = [
        v for v in variants.values()
        if not v.bam_injected and v.phase_set > 0 and
        v.hap1 in (0, 1) and v.hap2 in (0, 1) and v.hap1 != v.hap2
    ]
    left = max((v for v in anchors if v.pos <= gap_left), key=lambda v: v.pos)
    right = min((v for v in anchors if v.pos >= gap_right), key=lambda v: v.pos)
    return left.phase_set, right.phase_set


def truth_gauge(reads: dict[str, Read], truth: dict[str, str], phase_set: int):
    hap1_maternal = 0
    hap1_paternal = 0
    for read in reads.values():
        label = truth.get(read.name)
        if read.phase_set != phase_set or read.hap not in (1, 2) or label not in ("M", "P"):
            continue
        if (read.hap == 1 and label == "M") or (read.hap == 2 and label == "P"):
            hap1_maternal += 1
        else:
            hap1_paternal += 1
    if hap1_maternal == hap1_paternal:
        return None, hap1_maternal, hap1_paternal
    return hap1_maternal > hap1_paternal, hap1_maternal, hap1_paternal


def retained_component(
    variants: dict[int, Variant], observations: dict[str, dict[int, int]],
    left_ps: int, right_ps: int, all_categories: bool,
):
    allele_sets: dict[int, set[int]] = defaultdict(set)
    for row in observations.values():
        for site, allele in row.items():
            allele_sets[site].add(allele)

    eligible = {
        site for site, variant in variants.items()
        if (all_categories or (variant.category & GERMLINE_MASK) != 0) and
        allele_sets[site] == {0, 1}
    }
    left_sites = {
        site for site in eligible
        if variants[site].phase_set == left_ps and variants[site].hap1 != variants[site].hap2
    }
    right_sites = {
        site for site in eligible
        if variants[site].phase_set == right_ps and variants[site].hap1 != variants[site].hap2
    }

    parent = {site: site for site in eligible}

    def root(site: int) -> int:
        while parent[site] != site:
            parent[site] = parent[parent[site]]
            site = parent[site]
        return site

    def join(a: int, b: int):
        a = root(a)
        b = root(b)
        if a != b:
            parent[b] = a

    for row in observations.values():
        sites = [site for site in row if site in eligible]
        for site in sites[1:]:
            join(sites[0], site)

    roots_left = {root(site) for site in left_sites}
    roots_right = {root(site) for site in right_sites}
    bridge_roots = roots_left & roots_right
    retained = {site for site in eligible if root(site) in bridge_roots}
    return retained, left_sites & retained, right_sites & retained


def solve_orientation(
    variants: dict[int, Variant], observations: dict[str, dict[int, int]],
    sites: set[int], left_sites: set[int], right_sites: set[int], flip_right: bool,
    time_limit: float, weighted: bool,
):
    site_order = sorted(sites, key=lambda site: (variants[site].pos, site))
    site_col = {site: column for column, site in enumerate(site_order)}
    rows = {
        name: {site: allele for site, allele in row.items() if site in sites}
        for name, row in observations.items()
    }
    rows = {name: row for name, row in rows.items() if len(row) >= 2}
    read_order = sorted(rows)
    read_col = {
        name: len(site_order) + column for column, name in enumerate(read_order)
    }
    obs = [
        (name, site, allele)
        for name in read_order for site, allele in rows[name].items()
    ]
    error_offset = len(site_order) + len(read_order)
    variable_count = error_offset + len(obs)

    objective = np.zeros(variable_count)
    for obs_index, (_name, site, _allele) in enumerate(obs):
        objective[error_offset + obs_index] = (
            max(1, variants[site].weight) if weighted else 1
        )
    integrality = np.zeros(variable_count)
    integrality[:error_offset] = 1
    lower = np.zeros(variable_count)
    upper = np.ones(variable_count)

    for site in left_sites:
        orientation = variants[site].hap1
        lower[site_col[site]] = orientation
        upper[site_col[site]] = orientation
    for site in right_sites:
        orientation = variants[site].hap1 ^ int(flip_right)
        lower[site_col[site]] = orientation
        upper[site_col[site]] = orientation

    matrix_row = []
    matrix_col = []
    matrix_data = []
    constraint_lower = []
    constraint_upper = []

    def add_constraint(coefficients: list[tuple[int, float]], lb: float):
        row_index = len(constraint_lower)
        for column, value in coefficients:
            matrix_row.append(row_index)
            matrix_col.append(column)
            matrix_data.append(value)
        constraint_lower.append(lb)
        constraint_upper.append(np.inf)

    for obs_index, (name, site, allele) in enumerate(obs):
        x = site_col[site]
        y = read_col[name]
        error = error_offset + obs_index
        if allele == 0:
            add_constraint([(error, 1), (x, -1), (y, 1)], 0)
            add_constraint([(error, 1), (x, 1), (y, -1)], 0)
        else:
            add_constraint([(error, 1), (x, -1), (y, -1)], -1)
            add_constraint([(error, 1), (x, 1), (y, 1)], 1)

    matrix = coo_matrix(
        (matrix_data, (matrix_row, matrix_col)),
        shape=(len(constraint_lower), variable_count),
    ).tocsr()
    result = milp(
        c=objective,
        integrality=integrality,
        bounds=Bounds(lower, upper),
        constraints=LinearConstraint(matrix, constraint_lower, constraint_upper),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    assignments = {}
    read_costs = {}
    if result.x is not None:
        assignments = {name: int(result.x[column] >= 0.5) for name, column in read_col.items()}
        site_assignments = {site: int(result.x[column] >= 0.5) for site, column in site_col.items()}
        for name, row in rows.items():
            read_hap = assignments[name]
            read_costs[name] = sum(
                (max(1, variants[site].weight) if weighted else 1)
                for site, allele in row.items()
                if allele != (site_assignments[site] ^ read_hap)
            )
    return result, assignments, read_costs, len(site_order), len(read_order), len(obs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", type=Path)
    parser.add_argument("gap_left", type=int)
    parser.add_argument("gap_right", type=int)
    parser.add_argument("--truth", type=Path, default=Path("test_data/derived/chr20_truth_hap.tsv"))
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--weighted", action="store_true")
    parser.add_argument("--all-categories", action="store_true")
    parser.add_argument("--fold-mod", type=int, default=0)
    parser.add_argument("--fold-rem", type=int, default=0)
    args = parser.parse_args()

    variants, reads, observations = load_matrix(args.matrix)
    if args.fold_mod:
        observations = {
            name: row for name, row in observations.items()
            if int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little") % args.fold_mod == args.fold_rem
        }
    truth = {}
    with args.truth.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                truth[fields[0]] = fields[1][0]

    left_ps, right_ps = choose_flanks(variants, args.gap_left, args.gap_right)
    retained, left_sites, right_sites = retained_component(
        variants, observations, left_ps, right_ps, args.all_categories
    )
    left_gauge = truth_gauge(reads, truth, left_ps)
    right_gauge = truth_gauge(reads, truth, right_ps)
    truth_flip = None
    if left_gauge[0] is not None and right_gauge[0] is not None:
        truth_flip = left_gauge[0] != right_gauge[0]

    print(
        f"weighted={int(args.weighted)} all_categories={int(args.all_categories)} fold={args.fold_rem}/{args.fold_mod} gap={args.gap_left}-{args.gap_right} left_ps={left_ps} right_ps={right_ps} "
        f"sites={len(retained)} left_anchors={len(left_sites)} right_anchors={len(right_sites)} "
        f"truth_flip={truth_flip} gauges={left_gauge[1:]}/{right_gauge[1:]}"
    )
    if not retained or not left_sites or not right_sites:
        print("DISCONNECTED")
        return

    solved = []
    for flip in (False, True):
        result, assignments, read_costs, site_count, read_count, obs_count = solve_orientation(
            variants, observations, retained, left_sites, right_sites,
            flip, args.time_limit, args.weighted,
        )
        solved.append((result, assignments, read_costs))
        print(
            f"flip={int(flip)} status={result.status} success={result.success} "
            f"cost={result.fun} gap={getattr(result, 'mip_gap', None)} "
            f"sites={site_count} reads={read_count} observations={obs_count}"
        )

    if solved[0][0].fun is None or solved[1][0].fun is None:
        print("NO_SOLUTION")
        return
    if solved[0][0].fun == solved[1][0].fun:
        print("TIE")
        return
    predicted_flip = solved[1][0].fun < solved[0][0].fun
    winner = solved[int(predicted_flip)]
    shared_reads = solved[0][2].keys() & solved[1][2].keys()
    same_votes = sum(solved[0][2][name] < solved[1][2][name] for name in shared_reads)
    flip_votes = sum(solved[1][2][name] < solved[0][2][name] for name in shared_reads)
    vote_total = same_votes + flip_votes
    vote_p = binomtest(max(same_votes, flip_votes), vote_total, 0.5, alternative="greater").pvalue if vote_total else 1.0
    vote_flip = flip_votes > same_votes

    left_hap1_maternal = left_gauge[0]
    correct = total = 0
    if left_hap1_maternal is not None:
        for name, read_hap in winner[1].items():
            label = truth.get(name)
            if label not in ("M", "P"):
                continue
            predicted_maternal = (read_hap == 0) == left_hap1_maternal
            correct += predicted_maternal == (label == "M")
            total += 1
    print(
        f"WIN flip={int(predicted_flip)} delta={abs(solved[0][0].fun-solved[1][0].fun)} "
        f"votes={same_votes}/{flip_votes} vote_p={vote_p:.6g} vote_agree={vote_total > 0 and vote_flip == predicted_flip} "
        f"truth_match={truth_flip is None or predicted_flip == truth_flip} "
        f"read_truth={correct}/{total} ({correct/total if total else float('nan'):.4f})"
    )


if __name__ == "__main__":
    main()
