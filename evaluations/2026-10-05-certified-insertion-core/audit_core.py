#!/usr/bin/env python3
"""Require only the four diagnosed rescue reads to enter their existing core."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
expected = set(json.loads(Path('evaluations/2026-10-05-short-insertion-source-path/core-deficit-investigation.json').read_text()))

def read_output(folder):
    tags, votes = {}, defaultdict(Counter)
    with pysam.AlignmentFile(str(folder / 'phased.bam'), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            tags[read.query_name] = (hp, ps)
            if hp in (1, 2) and ps > 0 and read.query_name in truth:
                votes[ps][(hp == 1) != truth[read.query_name]] += 1
    orientation = {ps: c.most_common()[0][0] for ps, c in votes.items()}
    status = {q: ((hp == 1) != truth[q]) == orientation[ps]
              for q, (hp, ps) in tags.items() if hp in (1, 2) and ps > 0 and q in truth}
    return tags, status

before, before_status = read_output(args.before)
after, after_status = read_output(args.after)
changed = {q: {'before': before[q], 'after': after[q]}
           for q in before.keys() & after.keys() if before[q] != after[q]}
vcf_before = (args.before / 'phased.vcf').read_bytes()
vcf_after = (args.after / 'phased.vcf').read_bytes()
result = {'changed_tags': changed, 'same_output_molecules': before.keys() == after.keys(),
          'all_parental_status_preserved': before_status == after_status,
          'vcf_byte_identical': vcf_before == vcf_after,
          'vcf_sha256': hashlib.sha256(vcf_after).hexdigest()}
args.output.write_text(json.dumps(result, indent=2) + '\n')
assert result['same_output_molecules']
assert set(changed) == expected
for q in changed:
    old_hp, old_ps = before[q]
    hp, ps = after[q]
    assert old_hp == hp and old_ps == ps + 1_000_000_000 and ps > 0
    assert before_status[q] and after_status[q]
assert result['all_parental_status_preserved']
assert result['vcf_byte_identical']
print(json.dumps(result, indent=2))
