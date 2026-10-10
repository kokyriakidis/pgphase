#!/usr/bin/env python3
"""Compare the complete native/full HiPhase contract with the frozen step-5 audit."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
before = json.loads((ROOT / 'evaluations/2026-10-07-representation-step5/panel-audit.json').read_text())
after = json.loads((OUT / 'panel-audit.json').read_text())
old = {row['window']: row for row in before['windows']}
new = {row['window']: row for row in after['windows']}
assert old.keys() == new.keys()
changes = [dict(window=window, before=old[window], after=new[window])
           for window in old if old[window] != new[window]]
result = dict(panel_windows=len(new), changed_windows=len(changes),
              native_and_full_metrics_identical=not changes,
              summary_identical=before['summary'] == after['summary'], changes=changes)
(OUT / 'panel-comparison.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
