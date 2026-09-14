#!/usr/bin/env python3
"""Regression tests for full-sequence, multiallelic read observations."""
import unittest
import pysam
from allele_observations import observe_allele


def read(sequence, cigar, start=100):
    r = pysam.AlignedSegment()
    r.query_sequence = sequence
    r.reference_start = start
    r.cigartuples = cigar
    return r


class AlleleObservations(unittest.TestCase):
    def test_multiallelic_snp(self):
        self.assertEqual(observe_allele(read('AGT', [(0, 3)]), 101, 'C', ('C','G','T')), 1)

    def test_equal_length_insertions(self):
        alleles = ('A', 'ATT', 'AGG')
        self.assertEqual(observe_allele(read('AGGC', [(0,1),(1,2),(0,1)]),100,'A',alleles),2)
        self.assertIsNone(observe_allele(read('ATGC',[(0,1),(1,2),(0,1)]),100,'A',alleles))

    def test_multibase_substitution(self):
        self.assertIsNone(observe_allele(read('AC',[(0,2)]),100,'AA',('AA','AG','GT')))
        self.assertEqual(observe_allele(read('GT',[(0,2)]),100,'AA',('AA','AG','GT')),2)

    def test_deletion_and_partial_coverage(self):
        alleles=('ACT','A','AC')
        self.assertEqual(observe_allele(read('AG',[(0,1),(2,2),(0,1)]),100,'ACT',alleles),1)
        self.assertIsNone(observe_allele(read('AC',[(0,2)]),100,'ACT',alleles))

    def test_reference_skip_is_not_deletion(self):
        self.assertIsNone(observe_allele(read('AG',[(0,1),(3,2),(0,1)]),100,'ACT',('ACT','A')))

    def test_softclip_and_outside_insertion(self):
        self.assertEqual(observe_allele(read('TTAAGC',[(4,2),(0,1),(1,1),(0,2)]),101,'G',('G','T','C')),0)


if __name__ == '__main__':
    unittest.main()
