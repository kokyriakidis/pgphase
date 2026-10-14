#!/usr/bin/env python3
"""Every scored read group, ours against the rivals (evaluation only).

Usage: groups.py RUN_DIR

Reads: arms / centromere, each split by whether the read crosses a truth het
(phaseable) or not (nohet); each cell correct / discordant / unphased, as
score.py counts them.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402

GROUPS = ("arms_phaseable", "arms_nohet", "cen_phaseable", "cen_nohet", "all")


def main():
    runs = {"ours": sys.argv[1]}
    runs.update(D.competitors())
    print(f"{os.environ.get('PGDEV_CONTIG', '')}: correct / discordant / unphased")
    print(f"  {'tool':12}" + "".join(f"{g:>30}" for g in GROUPS))
    for name, d in runs.items():
        try:
            s = score.score_run(d)["reads"]
        except Exception as e:  # a rival without outputs for this profile
            print(f"  {name:12} {e}")
            continue
        print(f"  {name:12}" + "".join(
            f"{s[g]['correct']:>10,}/{s[g]['discordant']:>8,}/{s[g]['unphased']:>9,}"[-30:].rjust(30) for g in GROUPS))


if __name__ == "__main__":
    main()
