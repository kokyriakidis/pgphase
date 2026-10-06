#!/usr/bin/env python3
"""Measure HiPhase on the panel's identical primary alignments; cache measurements."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pysam
from cache_gap_replay import INDEX_SUFFIXES, identity


def alignment(read):
    return (read.reference_start, read.reference_end, read.cigarstring,
            hashlib.sha256(read.query_sequence.encode()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--hiphase', required=True)
    parser.add_argument('--truth', required=True)
    parser.add_argument('--panel', required=True)
    parser.add_argument('--out', default='src/test_gap_hiphase.tsv')
    args = parser.parse_args()
    output = Path(args.out)
    state_path = Path(args.out + '.json')
    inputs = []
    for path in (args.bam, args.hiphase, args.truth, args.panel):
        inputs.append(identity(path))
        for suffix in INDEX_SUFFIXES:
            index = Path(path + suffix)
            if index.exists():
                inputs.append(identity(index))
    signature = {'helper': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 'identity_helper': hashlib.sha256(
                     Path(__file__).with_name('cache_gap_replay.py').read_bytes()).hexdigest(),
                 'inputs': inputs}
    try:
        state = json.loads(state_path.read_text())
        if state['signature'] == signature and state['output_sha256'] == hashlib.sha256(output.read_bytes()).hexdigest():
            print('Reused HiPhase panel measurements')
            return
    except (OSError, ValueError, KeyError):
        pass
    truth = {fields[0]: fields[1].startswith('PATERNAL')
             for line in Path(args.truth).read_text().splitlines()
             if len(fields := line.split('\t')) == 2 and
             fields[1].startswith(('MATERNAL', 'PATERNAL'))}
    windows = []
    for line in Path(args.panel).read_text().splitlines():
        if not line or line.startswith('#') or line.startswith('gap_left'):
            continue
        fields = line.split('\t')
        windows.append((int(fields[0]), int(fields[1])))
    eligible, original = {}, {}
    with pysam.AlignmentFile(args.bam) as bam:
        for left, right in windows:
            names = set()
            for read in bam.fetch('CHM13#0#chr20', left - 1, right):
                if read.is_secondary or read.is_supplementary or read.query_name not in truth:
                    continue
                names.add(read.query_name)
                if read.query_name not in original:
                    original[read.query_name] = alignment(read)
            eligible[left, right] = names
    votes, tags = defaultdict(Counter), {}
    with pysam.AlignmentFile(args.hiphase) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            name = read.query_name
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            if name in truth and hp in (1, 2) and ps > 0:
                votes[ps][(hp == 1) != truth[name]] += 1
            if name in original:
                if alignment(read) != original[name]:
                    raise RuntimeError(f'HiPhase alignment differs from input: {name}')
                if name in tags:
                    raise RuntimeError(f'duplicate HiPhase primary read: {name}')
                tags[name] = hp, ps
    if tags.keys() != original.keys():
        raise RuntimeError('HiPhase does not contain every eligible input alignment')
    orientation = {ps: count.most_common(1)[0][0] for ps, count in votes.items()}
    lines = ['# Measured on identical primary input alignments, including unphased reads.',
             'window\tscorable\tcorrect\tcore_correct']
    for window, names in eligible.items():
        correct, core = 0, Counter()
        for name in names:
            hp, ps = tags[name]
            if hp in (1, 2) and ps in orientation and ((hp == 1) != truth[name]) == orientation[ps]:
                correct += 1
                core[ps] += 1
        lines.append(f'{window[0]}-{window[1]}\t{len(names)}\t{correct}\t{max(core.values(), default=0)}')
    output.write_text('\n'.join(lines) + '\n')
    state_path.write_text(json.dumps({'signature': signature,
                                    'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}, indent=2) + '\n')
    print(f'Measured {len(windows)} HiPhase windows on identical alignments')


if __name__ == '__main__':
    main()
