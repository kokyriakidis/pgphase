#!/usr/bin/env python3
"""Table for `pg sweep` (evaluation only).

stdin: lines "LABEL CHROM BASELINE_JSON RUN_JSON" (BASELINE_JSON "-" when the
profile has no pinned baseline) or "LABEL pg each: NAME failed". Per
configuration, each profile's change from its baseline and the total.
"""
import collections
import json
import sys

ROWS = [("correct", lambda s: s["reads"]["arms_phaseable"]["correct"]),
        ("discordant", lambda s: s["reads"]["arms_phaseable"]["discordant"]),
        ("unphased", lambda s: s["reads"]["arms_phaseable"]["unphased"]),
        ("switches", lambda s: s["variants"]["switches_arms"]),
        ("flips", lambda s: s["variants"]["flips_arms"]),
        ("hets phased", lambda s: s["variants"]["truth_matched_phased_hets"])]


def main():
    by_cfg = collections.OrderedDict()
    for line in sys.stdin:
        f = line.split()
        if not f:
            continue
        if "failed" in line:
            by_cfg.setdefault(f[0], []).append((" ".join(f[1:]), None, None))
            continue
        label, chrom, base, run = f
        by_cfg.setdefault(label, []).append((chrom, None if base == "-" else json.load(open(base)), json.load(open(run))))
    head = f"{'config':16} {'profile':8}" + "".join(f"{n:>22}" for n, _ in ROWS)
    print("arms_phaseable reads and arms variants; each cell: run (change from the pinned baseline)")
    print(head)
    for label, rows in by_cfg.items():
        tot = collections.Counter()
        have_base = True
        for chrom, base, run in rows:
            if run is None:
                print(f"{label:16} {chrom}")
                continue
            cells = []
            for name, get in ROWS:
                v = get(run)
                if base is None:
                    have_base = False
                    cells.append(f"{v:,}")
                else:
                    d = v - get(base)
                    tot[name] += d
                    cells.append(f"{v:,} ({d:+,})")
            print(f"{label:16} {chrom:8}" + "".join(f"{c:>22}" for c in cells))
        if have_base and len(rows) > 1:
            dc, dd = tot["correct"], tot["discordant"]
            ratio = f"   correct : discordant change {dc / dd:.1f} : 1" if dd > 0 and dc > 0 else ""
            print(f"{label:16} {'total':8}" + "".join(f"{tot[n]:>+22,}" for n, _ in ROWS) + ratio)


if __name__ == "__main__":
    main()
