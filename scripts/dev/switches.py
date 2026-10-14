#!/usr/bin/env python3
"""Every switch and flip on the arms, the run against the baseline (evaluation only).

Usage: switches.py RUN_DIR [BASELINE_DIR] [--all]

Exactly as score.py counts them: within each phase set, the truth-matched
phased hets in order; a change of orientation against the truth that reverts
at the next het is a flip, one that holds is a switch. Each event is printed
with its phase set, the last het before it and the first after it (a switch
lies between them), the hets on each side within the phase set, and whether
the baseline has an event between the same two hets ("both") or not ("new";
baseline-only events print as "gone"). --all prints flips too.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402


def events(run_dir):
    t = score.truth_tables()
    truth_phased = t["truth_phased"]
    vcf = os.path.join(run_dir, "phased.vcf")
    if not os.path.exists(vcf):
        vcf += ".gz"
    by_ps = collections.defaultdict(list)
    skip = score.excluded()
    for pos, ref, alts, gt, ps, cat in score.vcf_records(vcf):
        if ps is None or len(gt) != 2 or gt[0] == gt[1] or None in gt:
            continue
        if skip and score.in_regions(skip, pos - 1):
            continue
        if len(alts) == 1 and set(gt) == {0, 1}:
            key = (pos, ref.upper(), alts[0].upper())
            if key in truth_phased:
                by_ps[ps].append((pos, gt[0] == truth_phased[key]))
    out = []
    for ps, sites in by_ps.items():
        sites.sort()
        same = [x[1] for x in sites]
        i = 1
        while i < len(same):
            if same[i] != same[i - 1]:
                where = "cen" if score.CEN[0] <= sites[i][0] - 1 < score.CEN[1] else "arms"
                flip = i + 1 < len(same) and same[i + 1] == same[i - 1]
                out.append(dict(kind="flip" if flip else "switch", where=where, ps=ps, before=sites[i - 1][0],
                                after=sites[i][0], left=i, right=len(sites) - i,
                                block=(sites[0][0], sites[-1][0])))
                if flip:
                    i += 2
                    continue
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("baseline", nargs="?")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    kinds = ("switch", "flip") if a.all else ("switch",)
    run = [e for e in events(a.run) if e["where"] == "arms" and e["kind"] in kinds]
    base = [e for e in events(a.baseline) if e["where"] == "arms" and e["kind"] in kinds] if a.baseline else []
    bkeys = {(e["kind"], e["before"], e["after"]) for e in base}
    rkeys = {(e["kind"], e["before"], e["after"]) for e in run}
    rows = [(e, "both" if (e["kind"], e["before"], e["after"]) in bkeys else "new") for e in run]
    rows += [(e, "gone") for e in base if (e["kind"], e["before"], e["after"]) not in rkeys]
    rows.sort(key=lambda r: r[0]["before"])
    c = collections.Counter((e["kind"], s) for e, s in rows)
    print(f"{os.environ.get('PGDEV_CONTIG', '')}: arm " + ", ".join(
        f"{k}{'es' if k == 'switch' else 's'} {sum(v for (kk, s), v in c.items() if kk == k and s != 'gone')} "
        f"(new {c[(k, 'new')]}, gone {c[(k, 'gone')]})" for k in kinds))
    print(f"  {'status':6} {'kind':6} {'between hets':>27} {'hets left/right':>16} {'phase set span':>27}")
    for e, s in rows:
        print(f"  {s:6} {e['kind']:6} {e['before']:>13,}-{e['after']:<13,} {e['left']:>7}/{e['right']:<8} "
              f"{e['block'][0]:>13,}-{e['block'][1]:<13,}")


if __name__ == "__main__":
    main()
