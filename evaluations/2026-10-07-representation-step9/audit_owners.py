#!/usr/bin/env python3
"""Compare owning outputs and inspect retained graph/BAM joint-locus membership."""
from collections import Counter, defaultdict
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
WORK = ROOT / 'test_data/tmp_representation_step9/accepted'
spec = importlib.util.spec_from_file_location('read_audit',
    ROOT / 'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT / 'test_data/derived/chr20_truth_hap.tsv').read_text().splitlines()
         if len(f := line.split('\t')) == 2 and f[1] in ('PATERNAL', 'MATERNAL')}
results = []
for name in ['matrix4-final', 'matrix5-final', 'matrix65-final', 'whole-final']:
    before = ROOT / 'test_data/tmp_representation_step8/accepted' / name
    after = WORK / name
    bt, bs = audit.assignments(before / 'phased.bam', truth)
    at, ast = audit.assignments(after / 'phased.bam', truth)
    matrix = 'matrix.chunk0.bam-overlay-output.tsv'
    candidates_identical = (before / 'candidates.tsv').read_bytes() == (after / 'candidates.tsv').read_bytes()
    variants_identical = audit.variants(before / 'phased.vcf') == audit.variants(after / 'phased.vcf')
    matrix_identical = (before / matrix).read_bytes() == (after / matrix).read_bytes()
    for suffix in ['joint-contexts', 'joint-parents', 'joint-sequences', 'allele-contrasts', 'joint-molecules', 'joint-loci', 'joint-cohort', 'joint-alleles', 'joint-costs', 'joint-genotypes', 'joint-heldout']:
        assert (before / ('matrix.chunk0.' + suffix + '.tsv')).read_bytes() == (after / ('matrix.chunk0.' + suffix + '.tsv')).read_bytes(), suffix
    assert bt.keys() == at.keys() and bs.keys() == ast.keys()
    assert bt == at and bs == ast and candidates_identical and variants_identical and matrix_identical
    # Original observation rows must survive even when the clean solve changes.
    def observations(path):
        return [line for line in path.read_text().splitlines() if line.startswith('OBS\t')]
    assert observations(before / matrix) == observations(after / matrix)
    loci = defaultdict(list)
    for row in csv.DictReader((after / 'matrix.chunk0.joint-loci.tsv').open(), delimiter='\t'):
        loci[row['locus']].append(row)
    assert all(len(members) >= 2 and len({r['candidate'] for r in members}) == len(members)
               for members in loci.values())
    results.append(dict(before=str(before.relative_to(ROOT)), after=str(after.relative_to(ROOT)),
        candidate_rows_identical=candidates_identical, vcf_rows_identical=variants_identical,
        primary_tags_identical=bt == at, parental_status_identical=bs == ast,
        candidate_read_channel_quality_state_identical=matrix_identical,
        original_observations_identical=True, original_sequence_contexts_qualities_and_provenance_identical=True,
        changed_primary_tags=sum(bt[name] != at[name] for name in bt),
        parental_transitions=dict(Counter(f'{bs[name]}->{ast[name]}' for name in bs if bs[name] != ast[name])),
        before_parental_counts=dict(Counter(bs.values())),
        matrix_sha256=hashlib.sha256((after / matrix).read_bytes()).hexdigest(),
        primary_reads=len(at), parental_counts=dict(Counter(ast.values())),
        equivalent_loci=len(loci), retained_alias_rows=sum(map(len, loci.values())),
        mixed_graph_bam_loci=sum({r['bam_injected'] for r in members} == {'0', '1'}
                                for members in loci.values()),
        mixed_categories=sum(len({r['category'] for r in members}) > 1 for members in loci.values()),
        loci=list(loci.values())))
(OUT / 'owner-checks.json').write_text(json.dumps(results, indent=2) + '\n')
for report in results:
    print(json.dumps({k: v for k, v in report.items() if k != 'loci'}, indent=2))
