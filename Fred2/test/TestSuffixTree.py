import numpy as np
from unittest import TestCase

from Fred2.Utility import generate_overlap_graph


class TestOverlappingGrapg(TestCase):

    def test_overlapping_graph(self):
        sequences = ['ACGCA', 'CGCAT',"TCGCG", "CATTC","ATTCG"]
        adj = generate_overlap_graph(sequences)
        expected = np.array([[0, 5, 5, 5, 5, 5],
                             [0, 0, 1, 5, 3, 4],
                             [0, 5, 0, 4, 2, 3],
                             [0, 5, 3, 0, 5, 5],
                             [0, 5, 4, 3, 0, 1],
                             [0, 5, 3, 2, 5, 0]])
        self.assertTrue((adj==expected).all())




