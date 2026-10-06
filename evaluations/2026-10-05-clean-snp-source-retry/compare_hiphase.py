#!/usr/bin/env python3
"""Compare the same original reads in the gap, including abstentions."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
left, right = 57764235, 57785224
eligible = {}
input_names = set()
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
    for read in bam.fetch('CHM13#0#chr20', left - 1, right):
        if read.is_secondary or read.is_supplementary:
            continue
        input_names.add(read.query_name)
        if read.query_name not in truth:
            continue
        eligible[read.query_name] = (read.reference_start, read.reference_end,
                                     read.cigarstring, read.query_sequence)
assert input_names == eligible.keys()
outputs = [('pgphase_before', 'test_data/tmp_gap_fix63/full-deferred-source-runs/phased.bam'),
           ('pgphase_after', 'test_data/tmp_gap_fix64/full-source-retry/phased.bam'),
           ('hiphase', 'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam')]
reports, statuses = [], {}
for name, path in outputs:
    votes, tags, same_alignments = defaultdict(Counter), {}, 0
    with pysam.AlignmentFile(path, check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            if hp in (1, 2) and ps > 0 and qname in truth:
                votes[ps][(hp == 1) != truth[qname]] += 1
            if qname in eligible:
                tags[qname] = (hp, ps)
                if name == 'hiphase':
                    same_alignments += eligible[qname] == (read.reference_start, read.reference_end,
                                                          read.cigarstring, read.query_sequence)
    orientation = {ps: c.most_common()[0][0] for ps, c in votes.items()}
    status, local = {}, defaultdict(Counter)
    for qname in eligible:
        hp, ps = tags.get(qname, (0, 0))
        if hp not in (1, 2) or ps not in orientation:
            status[qname] = 'unphased'
        else:
            correct = ((hp == 1) != truth[qname]) == orientation[ps]
            status[qname] = 'correct' if correct else 'discordant'
            local[ps][status[qname]] += 1
    statuses[name] = status
    counts = Counter(status.values())
    reports.append({'stage': name, 'source': path, 'input_primary_reads': len(input_names), 'truth_scorable': len(eligible),
                    'counts': dict(counts), 'correct_all_reads': counts['correct'] / len(eligible),
                    'read_phase_sets': len(local),
                    'dominant_correct': max((c['correct'] for ps, c in local.items()
                                             if name == 'hiphase' or ps < 1_000_000_000), default=0),
                    'identical_input_alignments': same_alignments if name == 'hiphase' else None})
    if name == 'hiphase':
        assert same_alignments == len(eligible)
result = {'measurements': reports,
          'hiphase_correct_pgphase_not': [q for q in eligible if statuses['hiphase'][q] == 'correct' and statuses['pgphase_after'][q] != 'correct'],
          'pgphase_correct_hiphase_not': [q for q in eligible if statuses['pgphase_after'][q] == 'correct' and statuses['hiphase'][q] != 'correct']}
(Path('evaluations/2026-10-05-clean-snp-source-retry') / 'hiphase-comparison.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(reports, indent=2))

pg, hp = reports[1], reports[2]
assert pg["correct_all_reads"] >= 0.80
assert pg["counts"]["correct"] >= hp["counts"]["correct"]
assert pg["dominant_correct"] >= hp["dominant_correct"]
assert pg["counts"].get("discordant", 0) <= hp["counts"].get("discordant", 0)
