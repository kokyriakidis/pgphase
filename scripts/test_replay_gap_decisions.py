#!/usr/bin/env python3
import tempfile
import json
from pathlib import Path
import unittest

from replay_gap_decisions import classify, orientation_costs, replay, summarize, resolve_components, direct_observation_support, prepare_gap, score_gap, load_features, haplotype_margin, FEATURE_SCHEMA


class ReplayTests(unittest.TestCase):
    def test_parity_composition_is_order_independent(self):
        edges = [(100, 200, True), (200, 300, True), (100, 300, False)]
        a = resolve_components(edges)
        self.assertEqual(a, resolve_components(list(reversed(edges))))
        self.assertTrue(a[0]['consistent'])
        self.assertEqual(a[0]['flips'], {'100': False, '200': True, '300': False})

    def test_conflicting_component_has_no_assignments(self):
        result = resolve_components([(100, 200, True), (200, 300, True), (100, 300, True)])
        self.assertFalse(result[0]['consistent'])
        self.assertEqual(result[0]['flips'], {})
        self.assertEqual(result[0]['nodes'], [100, 200, 300])

    def test_untagged_bridge_and_correlated_events(self):
        gap = [0, 100, 300, 150, 250, 50, 350]
        sites = {0: dict(pos=100, type=8, length=1, ps=100, h1=0, h2=1),
                 1: dict(pos=300, type=8, length=1, ps=300, h1=0, h2=1),
                 2: dict(pos=100, type=2, length=1, ps=100, h1=0, h2=1),
                 3: dict(pos=99, type=8, length=1, ps=100, h1=0, h2=1)}
        reads = {0: dict(molecule='0:bridge', hp=0, skipped=False, flag=0)}
        observations = [(0, 3, 0, 0, 0, 40), (0, 0, 0, 0, 0, 40),
                        (0, 1, 1, 1, 1, 40), (0, 2, 0, 0, 0, 40)]
        result = direct_observation_support(gap, sites, reads, observations)
        self.assertEqual(result['hybrid']['flip'], 1)
        self.assertEqual(result['hybrid']['bridges'][0]['counts'], [[1, 0], [0, 1]])
        # Conflicting descriptions of one local event are not independent votes.
        observations[-1] = (0, 2, 1, 0, 0, 40)
        result = direct_observation_support(gap, sites, reads, observations)
        self.assertEqual(result['hybrid']['flip'], 0)
        self.assertEqual(result['graph']['flip'], 1)

    def test_large_total_does_not_hide_missing_haplotype_support(self):
        bad = [[1, 0, 0, 2], [0, 15, 17, 0]]
        self.assertEqual(classify(bad, 2), 'FLIP')
        self.assertEqual(haplotype_margin(bad, True), 1)
        good = [[0, 22, 26, 0], [0, 32, 26, 0]]
        self.assertEqual(haplotype_margin(good, False), 22)

    def test_missing_is_not_opposition(self):
        self.assertEqual(classify([[20, 0, 0, 20], [0, 0, 0, 0]], 2), 'INSUFFICIENT')
        self.assertEqual(classify([[5, 5, 5, 5], [5, 5, 5, 5]], 2), 'CONFLICTING')

    def test_label_gauge(self):
        votes = [[20, 1, 2, 10], [1, 15, 12, 2]]
        self.assertEqual(classify(votes, 2), 'FLIP')
        swapped = [[v[1], v[0], v[3], v[2]] for v in votes]
        self.assertEqual(orientation_costs(votes), orientation_costs(swapped))
        old_right_swapped = [votes[0], votes[1][2:] + votes[1][:2]]
        self.assertEqual(classify(old_right_swapped, 2), 'SAME')

    def test_influence_distinguishes_abstention(self):
        result = summarize([[2, 0, 0, 0], [20, 0, 0, 0]],
                           [('weak-flank-read', 0, 0, 1)], 2)
        self.assertEqual(result['decision'], 'SAME')
        self.assertEqual(result['support_critical_molecules'], ['weak-flank-read'])
        self.assertEqual(result['reversing_molecules'], [])

    def fixture(self, root):
        prefix = root / 'gap'
        Path(str(prefix) + '.evidence.tsv').write_text(
            'SCHEMA\t1\nGAP\t0\t100\t300\t100\t300\t50\t350\n'
            'READ\t0\t0\tleft\t50\t200\t60\t0\t1\t100\t0\n'
            'READ\t1\t0\tright\t200\t350\t60\t0\t1\t300\t0\n')
        Path(str(prefix) + '.reads.tsv').write_text(
            'VIEW\tTIER\tREAD\tHP\tPS\tSKIPPED\n'
            '0\t1\t0\t1\t200\t0\n0\t1\t1\t1\t200\t0\n'
            '1\t1\t0\t1\t200\t0\n1\t1\t1\t2\t200\t0\n')
        Path(str(prefix) + '.votes.tsv').write_text(
            'VIEW\tTIER\tPS\tL11\tL12\tL21\tL22\tR11\tR12\tR21\tR22\tJOINED\tFLIP\n'
            '0\t1\t200\t1\t0\t0\t0\t1\t0\t0\t0\t0\t0\n'
            '1\t1\t200\t1\t0\t0\t0\t0\t1\t0\t0\t0\t0\n')
        return prefix

    def test_views_are_alternatives_not_added_votes(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = self.fixture(Path(tmp))
            self.assertEqual(replay(prefix, margin=1)['decision'], 'CONFLICTING')
            self.assertEqual(replay(prefix, margin=2)['decision'], 'INSUFFICIENT')

    def test_replay_detects_corrupted_vote_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = self.fixture(Path(tmp))
            path = Path(str(prefix) + '.reads.tsv')
            path.write_text(path.read_text().replace('0\t1\t0\t1\t200\t0', '0\t1\t0\t2\t200\t0'))
            with self.assertRaisesRegex(ValueError, 'Exported votes differ'):
                replay(prefix)

    def test_feature_cache_reuses_inputs_and_rejects_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = self.fixture(root)
            prepared = prepare_gap(prefix)
            cache = root / 'features.json'
            cache.write_text(json.dumps(dict(schema=FEATURE_SCHEMA, gaps=[prepared])))
            paths = [Path(str(prefix) + '.evidence.tsv')]
            loaded = load_features(cache, paths)[0]
            self.assertEqual(score_gap(loaded, 1), replay(prefix, 1))
            self.assertEqual(score_gap(loaded, 2), replay(prefix, 2))
            paths[0].write_text(paths[0].read_text() + '\n')
            with self.assertRaisesRegex(ValueError, 'Frozen input changed'):
                load_features(cache, paths)

    def test_duplicate_molecules_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = self.fixture(Path(tmp))
            path = Path(str(prefix) + '.evidence.tsv')
            path.write_text(path.read_text().replace('\tright\t', '\tleft\t'))
            with self.assertRaisesRegex(ValueError, 'Duplicate molecule'):
                replay(prefix)


if __name__ == '__main__':
    unittest.main()
