#!/usr/bin/env python3
"""Compare full-chromosome output and parental read assignments with the starting state."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BEFORE = ROOT / "test_data/tmp_representation_step5/accepted/full-current"
AFTER = ROOT / "test_data/tmp_representation_step6/accepted/full-current"
spec = importlib.util.spec_from_file_location("read_audit",
    ROOT / "evaluations/2026-10-06-shared-insertion-source/read_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {fields[0]: fields[1] == "PATERNAL"
         for line in (ROOT / "test_data/derived/chr20_truth_hap.tsv").read_text().splitlines()
         if len(fields := line.split("\t")) == 2 and
         fields[1] in ("PATERNAL", "MATERNAL")}
old_variants = audit.variants(BEFORE / "phased.vcf")
new_variants = audit.variants(AFTER / "phased.vcf")
old_tags, old_status = audit.assignments(BEFORE / "phased.bam", truth)
new_tags, new_status = audit.assignments(AFTER / "phased.bam", truth)
assert old_tags.keys() == new_tags.keys()
assert old_status.keys() == new_status.keys()
changed_tags = [name for name in old_tags if old_tags[name] != new_tags[name]]
changed_status = [name for name in old_status if old_status[name] != new_status[name]]
transitions = Counter(f'{old_status[name]}->{new_status[name]}' for name in changed_status)
previously_phased_unchanged = all(old_tags[name] == new_tags[name] for name in old_tags
                                if old_tags[name][0] in (1, 2) and old_tags[name][1] > 0)
core_unchanged = all(old_tags[name] == new_tags[name] for name in old_tags
                     if old_tags[name][0] in (1, 2) and 0 < old_tags[name][1] < 1_000_000_000)


def blocks(variants):
    lengths = sorted((hi - lo + 1 for lo, hi in audit.block_extents(variants)), reverse=True)
    accumulated = 0
    n50 = 0
    for length in lengths:
        accumulated += length
        if accumulated >= sum(lengths) / 2:
            n50 = length
            break
    return dict(count=len(lengths), n50_bp=n50, largest_block_bp=max(lengths, default=0))


result = dict(binary_sha256=hashlib.sha256((ROOT / "pgphase").read_bytes()).hexdigest(),
              candidates_byte_identical=(BEFORE / "candidates.tsv").read_bytes() ==
                                        (AFTER / "candidates.tsv").read_bytes(),
              vcf_rows_identical=old_variants == new_variants,
              before_vcf_rows=len(old_variants), after_vcf_rows=len(new_variants),
              primary_reads=len(old_tags), scorable_primary_reads=len(old_status),
              changed_primary_hp_ps=len(changed_tags),
              changed_parental_status=len(changed_status), changed_status_transitions=dict(transitions), previously_phased_tags_unchanged=previously_phased_unchanged, connected_core_tags_unchanged=core_unchanged,
              before_parental_counts=dict(Counter(old_status.values())),
              after_parental_counts=dict(Counter(new_status.values())),
              before_blocks=blocks(old_variants), after_blocks=blocks(new_variants),
              changed_primary_examples=[dict(read=name, before=old_tags[name], after=new_tags[name])
                                        for name in sorted(changed_tags)])
print(json.dumps(result, indent=2))
