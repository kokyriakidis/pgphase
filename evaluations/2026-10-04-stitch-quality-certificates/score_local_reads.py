"""Score unique primary molecules physically overlapping each gap after phasing."""
import json
from pathlib import Path
from collections import Counter, defaultdict
import pysam
root = Path('test_data/tmp_gap_fix50')
truth = {f[0]: f[1] for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
results = []
for mb, left, right in [(35, 35498368, 35516845), (57, 57854341, 57866713)]:
    with pysam.AlignmentFile(str(root / 'hiphase' / str(mb) / 'input.bam')) as bam:
        eligible = {r.query_name for r in bam
                    if not (r.is_secondary or r.is_supplementary) and
                    r.reference_start <= right and r.reference_end > left}
    for stage, path in [('pgphase', root / 'baseline' / str(mb) / 'phased.bam'),
                        ('hiphase', root / 'hiphase' / str(mb) / 'phased.bam')]:
        votes = defaultdict(Counter)
        seen = set()
        with pysam.AlignmentFile(str(path), check_sq=False) as bam:
            for read in bam:
                name = read.query_name
                if (read.is_secondary or read.is_supplementary or
                    name in seen or name not in eligible or name not in truth):
                    continue
                seen.add(name)
                hp = read.get_tag('HP') if read.has_tag('HP') else 0
                ps = read.get_tag('PS') if read.has_tag('PS') else 0
                if hp in (1, 2) and ps > 0:
                    votes[ps][(hp == 1) == (truth[name] == 'PATERNAL')] += 1
        scored = sum(sum(c.values()) for c in votes.values())
        correct = sum(max(c.values()) for c in votes.values())
        results.append(dict(chunk=mb, stage=stage, gap=[left, right], scored=scored,
                            correct=correct, discordant=scored - correct,
                            phase_sets=len(votes), accuracy=correct / scored if scored else None))
output = Path('evaluations/2026-10-04-stitch-quality-certificates/local-read-truth.json')
output.write_text(json.dumps(results, indent=2) + '\n')
print(output.read_text())
