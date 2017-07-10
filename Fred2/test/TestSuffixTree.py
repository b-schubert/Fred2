import numpy as np
import multiprocessing as mp
from unittest import TestCase

from Fred2.Core import Peptide
from Fred2.EpitopeSelection.Mosaic import _suffixPrefixMatch
from Fred2.Utility import GeneralisedSuffixTree, SuffixTree, generate_overlap_graph


class TestPeptide(TestCase):

    # def setUp(self):
    #     self.simple = Peptide("callage")
    #     self.simple = Peptide("mississippi")
    #     self.simpleSet = [Peptide("callage"),Peptide("aac")]
    #     self.combined = Peptide("callageaac")



    def _st(self):
        sequences = ['ACGCA', 'CGCAT', "TCGCG", "CATTC", "ATTCG"]
        st = GeneralisedSuffixTree(sequences)
        for n in st.preOrderNodes:
            p = n.parent
            p_rank = -1 if p is None else p.dfsRank
            p_idx = -1 if p is None else p.index
            print("P_index",p_rank,"EdgeLabel",n.edgeLabel,"Index",n.dfsRank,
                  "TreeDepth",n.treeDepth, "DFSRank", n.dfsRank, "isLeaf", n.isLeaf)

    def test_overlapping_graph(self):
        sequences = ['ACGCA', 'CGCAT',"TCGCG", "CATTC","ATTCG"]
        adj = generate_overlap_graph(sequences)
        expected = np.array([[0, 4, 0, 2, 1],
                             [0, 0, 1, 3, 2],
                             [0, 2, 0, 0, 0],
                             [0, 1, 2, 0, 4],
                             [0, 2, 3, 0, 0]])
        self.assertTrue((adj==expected).all())

    # def test_simple_suffix_tree(self):
    #     print 'SIMPLE TEST'
    #     st = SuffixTree('mississippi', '#')
    #     assert st.string == 'mississippi#'
    #     st = SuffixTree('mississippi')
    #     assert st.string == 'mississippi$'
    #
    #     r = st.root
    #     assert st.root == r
    #     assert st.root.parent is None
    #     assert st.root.firstChild.parent is not None
    #     assert st.root.firstChild.parent == st.root
    #
    #     for n in st.postOrderNodes:
    #         assert st.string[n.start:n.end + 1] == n.edgeLabel
    #
    #     # collect path labels
    #     for n in st.preOrderNodes:
    #         p = n.parent
    #         if p is None:  # the root
    #             n._pathLabel = ''
    #         else:
    #             n._pathLabel = p._pathLabel + n.edgeLabel
    #
    #     for n in st.postOrderNodes:
    #         assert n.pathLabel == n._pathLabel
    #
    #     for l in st.leaves:
    #         print 'leaf:', '"' + l.pathLabel + '"', ':', '"' + l.edgeLabel + '"'
    #
    #     for n in st.innerNodes:
    #         print 'inner:', '"' + n.edgeLabel + '"'
    #     print 'done.\n\n'
    #
    #     del st
    #
    # def test_generalized_suffix_tree(self):
    #     print 'GENERALISED TEST'
    #     sequences = ['xabxa', 'babxba']
    #     st = GeneralisedSuffixTree(sequences)
    #
    #     for shared in st.sharedSubstrings():
    #         print '-' * 70
    #         for seq, start, stop in shared:
    #             print seq, '[' + str(start) + ':' + str(stop) + ']',
    #             print sequences[seq][start:stop],
    #             print sequences[seq][:start] + '|' + sequences[seq][start:stop] + \
    #                   '|' + sequences[seq][stop:]
    #     print '=' * 70
    #
    #     for shared in st.sharedSubstrings(2):
    #         print '-' * 70
    #         for seq, start, stop in shared:
    #             print seq, '[' + str(start) + ':' + str(stop) + ']',
    #             print sequences[seq][start:stop],
    #             print sequences[seq][:start] + '|' + sequences[seq][start:stop] + \
    #                   '|' + sequences[seq][stop:]
    #     print '=' * 70
    #
    #     print 'done.\n\n'



