#!/usr/bin/env python3
"""Score controlled HiPhase inputs and audit named read segments at gap sites."""
import argparse
import hashlib
import importlib.util
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--replays', type=Path, default=ROOT/'test_data/tmp_joint_evidence_investigation')
args = parser.parse_args()
directory = args.replays/'hiphase5'
spec = importlib.util.spec_from_file_location('read_audit', ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
with pysam.AlignmentFile(str(directory/'input.bam')) as bam:
    names = {r.query_name for r in bam.fetch('CHM13#0#chr20', 5309405, 5345085)
             if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
def geometry(path):
    with pysam.AlignmentFile(str(path)) as bam:
        return {r.query_name: (r.reference_name, r.reference_start, r.reference_end, r.cigarstring,
                               r.flag, r.mapping_quality, r.query_sequence, bytes(r.query_qualities))
                for r in bam if not r.is_secondary and not r.is_supplementary}
input_geometry = geometry(directory/'input.bam')
pg_variants, pg_calls = {}, defaultdict(dict)
for line in (args.replays/'5/matrix.chunk0.bam-overlay-output.tsv').open():
    fields = line.rstrip().split('\t')
    if fields[0] == 'VAR':
        pg_variants[int(fields[1])] = fields
    elif fields[0] == 'OBS':
        pg_calls[fields[1]][int(fields[2])] = int(fields[3])
pg_indices = {pos: next(i for i, v in pg_variants.items() if int(v[2]) == pos and v[3] == variant_type)
              for pos, variant_type in [(5309406, 'X'), (5315591, 'I'), (5331266, 'D'), (5345085, 'X')]}
pg_pairs = {}
for left, right in [(5309406, 5315591), (5315591, 5331266)]:
    pg_pairs[f'{left}-{right}'] = sum(p.get(pg_indices[left]) in (0, 1) and p.get(pg_indices[right]) in (0, 1) for q, p in pg_calls.items() if q in names)
report = {'gap': [5309406, 5345085], 'scorable': len(names), 'pgphase_callable_pairs': pg_pairs, 'controls': {}}
for mode in ['default', 'local', 'dv', 'add-marker', 'diploid', 'both']:
    input_path = directory/('input.vcf.gz' if mode in ('default', 'local') else mode+'-input.vcf.gz')
    with pysam.VariantFile(str(input_path)) as vcf:
        variants = [r.copy() for r in vcf if (gt := next(iter(r.samples.values())).get('GT')) and None not in gt and len(set(gt)) == 2]
    blocks = []
    for line in (directory/(mode+'.log')).open():
        if 'Solving problem: PhaseBlock' in line:
            match = re.search(r'coordinates: "(.+):(\d+)-(\d+)".*num_variants: (\d+)', line)
            assert match, line
            blocks.append({'beg': int(match[2]), 'end': int(match[3]), 'nvars': int(match[4]), 'rows': []})
        elif 'read segment #' in line and 'read_name:' in line:
            match = re.search(r'read_name: "([^"]+)", alleles: \[(.*?)\], quals: \[(.*?)\], region: (\d+)\.\.(\d+)', line)
            assert match, line
            alleles = match[2].split(', ')
            beg, end = int(match[4]), int(match[5])
            assert len(alleles) == end-beg
            blocks[-1]['rows'].append((match[1], beg, end, alleles))
    calls = defaultdict(dict)
    for block in blocks:
        sites = [r for r in variants if block['beg'] <= r.start <= block['end']]
        assert len(sites) == block['nvars'], (mode, len(sites), block['nvars'])
        for q, beg, end, alleles in block['rows']:
            assert end <= len(sites)
            for i, allele in enumerate(alleles, beg):
                if allele in ('Reference', 'Alternate'):
                    calls[q][sites[i].pos] = 0 if allele == 'Reference' else 1
    tags, status = audit.assignments(directory/(mode+'.bam'), truth)
    assert geometry(directory/(mode+'.bam')) == input_geometry, mode
    correct_cores = Counter(tags[q][1] for q in names if status.get(q) == 'correct' and 0 < tags[q][1] < 1_000_000_000)
    pairs = {f'{left}-{right}': sum(left in calls[q] and right in calls[q] for q in names)
             for left, right in [(5309406, 5315591), (5315591, 5331265), (5331265, 5339363), (5339363, 5345085)]}
    elapsed = re.search(r'after ([\d.]+) seconds', (directory/(mode+'.log')).read_text())
    report['controls'][mode] = {'input_sha256': hashlib.sha256(input_path.read_bytes()).hexdigest(),
                                'counts': dict(Counter(status.get(q, 'unphased') for q in names)),
                                'core_correct': max(correct_cores.values(), default=0), 'correct_cores': dict(correct_cores),
                                'phased_blocks': len(blocks), 'callable_pairs': pairs,
                                'identical_primary_alignments': len(input_geometry),
                                'seconds': float(elapsed[1]) if elapsed else None}
    print(mode, report['controls'][mode])
(OUT/'competitor-controls.json').write_text(json.dumps(report, indent=2)+'\n')
