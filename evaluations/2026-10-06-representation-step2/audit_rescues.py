#!/usr/bin/env python3
"""Measure the newly normalized output rescues without treating them as connected core."""
from collections import Counter
import importlib.util
import json
from pathlib import Path

import pysam
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('read_audit', ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').read_text().splitlines() if len(f := line.split('\t')) == 2 and f[1] in ('PATERNAL', 'MATERNAL')}
old_tags, old_status = audit.assignments(ROOT/'test_data/tmp_representation_step1/full-current/phased.bam', truth)
new_tags, new_status = audit.assignments(ROOT/'test_data/tmp_representation_step2/full-current/phased.bam', truth)
hi_tags, hi_status = json.loads((ROOT/'test_data/tmp_gap_fix78/hi-tags.json').read_text())
changed = sorted(q for q in old_tags if old_tags[q] != new_tags[q])
rows = {q: dict(read=q, before=old_tags[q], after=new_tags[q], parental_status=new_status[q], hiphase=hi_tags.get(q), hiphase_parental_status=hi_status.get(q, 'unphased')) for q in changed}
if rows:
    with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
        for read in bam:
            if read.query_name in rows and not read.is_secondary and not read.is_supplementary:
                rows[read.query_name].update(beg=read.reference_start + 1, end=read.reference_end, mapq=read.mapping_quality)

assert all('beg' in row for row in rows.values())
report = dict(changed_primary_reads=len(rows), pgphase_counts=dict(Counter(new_status[q] for q in changed)), hiphase_counts=dict(Counter(hi_status.get(q, 'unphased') for q in changed)), output_rescue_only=all(new_tags[q][1] >= 1_000_000_000 for q in changed), status_pairs=dict(Counter(f'{new_status[q]}/{hi_status.get(q, "unphased")}' for q in changed)), reads=list(rows.values()))
(OUT/'rescues-checks.json').write_text(json.dumps(report, indent=2)+'\n')
print({k: v for k, v in report.items() if k != 'reads'})
print([row for row in rows.values() if row['parental_status'] == 'discordant'])
