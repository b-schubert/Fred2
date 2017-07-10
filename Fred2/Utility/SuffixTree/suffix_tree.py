import re
import numpy as np
import Fred2._suffix_tree


def postOrderNodes(node):
    '''Iterator through all nodes in the sub-tree rooted in node in
    post-order.'''
    def dfs(n):
        c = n.firstChild
        while c is not None:
            for m in dfs(c):
                yield m
            c = c.next
        yield n
    for n in dfs(node):
        yield n


def preOrderNodes(node):
    '''Iterator through all nodes in the sub-tree rooted in node in
    pre-order.'''
    def dfs(n):
        yield n
        c = n.firstChild
        while c is not None:
            for m in dfs(c):
                yield m
            c = c.next
    for n in dfs(node):
        yield n


def leaves(node):
    'Iterator through all leaves in the tree rooted in node.'
    for n in postOrderNodes(node):
        if n.isLeaf:
            yield n


def innerNodes(node):
    'Iterator through all inner nodes in the tree rooted in node.'
    for n in postOrderNodes(node):
        if not n.isLeaf:
            yield n


def children(node):
    'Iterate through all immediate children of node.'
    c = node.firstChild
    while c is not None:
        yield c
        c = c.next


class SuffixTree(Fred2._suffix_tree.SuffixTree):

    """A higher-level wrapper around the C suffix tree type,
_suffix_tree.SuffixTree.  This class adds a few methods to the suffix
tree, methods that are more easily expressed in Python than in C, and
that can be written using the primitives exported from C.  """

    def __init__(self,s,t='$'):
        '''Build a suffix tree from the input string s. The string
must not contain the special symbol $.'''
        if t in s:
            raise "The suffix tree string must not contain terminal symbol!"
        Fred2._suffix_tree.SuffixTree.__init__(self,s,t)
        self.label_tree_depth()
        for i,n in enumerate(self.preOrderNodes):
            n.dfsRank = i

    def label_tree_depth(self):
        'Label deapth of the tree TODO: MOVE TO C CODE'
        def dfs(n, depth):
            n.treeDepth = depth
            c = n.firstChild
            while c is not None:
                dfs(c, depth+1)
                c = c.next
        dfs(self.root, 0)

    def generatePostOrderNodes(self):
        'Iterator through all nodes in the tree in post-order.'
        for n in postOrderNodes(self.root):
            yield n

    def generatePreOrderNodes(self):
        'Iterator through all nodes in the tree in pre-order.'
        for n in preOrderNodes(self.root):
            yield n

    def generateLeaves(self):
        'Iterator through all leaves in the tree.'
        for n in leaves(self.root):
            yield n

    def generateInnerNodes(self):
        'Iterator through all leaves in the tree.'
        for n in innerNodes(self.root):
            yield n

    # set class properties
    postOrderNodes = property(generatePostOrderNodes, None, None,
                              "postOrderNodes")
    preOrderNodes = property(generatePreOrderNodes, None, None,
                             "preOrderNodes")
    
    leaves = property(generateLeaves, None, None, "leaves")
    innerNodes = property(generateInnerNodes, None, None, "innerNodes")


class GeneralisedSuffixTree(SuffixTree):

    """A suffix tree for a set of strings."""

    def __init__(self, sequences):
        '''Build a generalised suffix tree.  The strings must not
contain the special symbols $ or ascii numbers from 1 to the number of
sequences.'''

        self.sequences = sequences
        self.startPositions = [0]
        concatString = ''
        for i in xrange(len(sequences)):
            if chr(i+1) in sequences[i]:
                raise ValueError("The suffix tree string must not contain chr(%d)!"%(i+1))
            concatString += str(sequences[i])+str(i+1)
            self.startPositions += [len(concatString)]

        self.startPositions += [self.startPositions[-1]+1] # empty string
        self.sequences += ['']

        SuffixTree.__init__(self, concatString)
        self.__annotateNodes()

    def _translateIndex(self, idx):
        'Translate a concat-string index into a (stringNo,idx) pair.'
        for i in xrange(len(self.startPositions)-1):
            if self.startPositions[i] <= idx < self.startPositions[i+1]:
                return (i, idx-self.startPositions[i])
        raise IndexError("Index out of range: "+str(idx))

    def __annotateNodes(self):
        for n in self.postOrderNodes:
            if n.isLeaf:
                seq,idx = self._translateIndex(n.index)
                n.pathIndices = [(seq, idx)]
                n.sequences = [seq]
                n.L = []

            else:
                pathIndices = [] ; sequences = []; L = []
                c = n.firstChild
                while c is not None:
                    pathIndices += c.pathIndices
                    sequences += c.sequences
                    idx = c.edgeLabel[0]
                    print(n.dfsRank, c.isLeaf, idx, set(sequences))
                    if c.isLeaf and idx.isdigit() and n.start != -1:
                        L.append(int(idx)-1)
                    c = c.next

                seqsInSubtree = {}
                for s in sequences:
                    seqsInSubtree[s] = 1

                n.pathIndices = pathIndices
                n.sequences = [s for s in seqsInSubtree]
                n.L = list(set(L))

    def sharedSubstrings(self, minimumLength=0):
        '''Iterator through shared sub-strings.  Returned as a list of triples
 (sequence,from,to).'''
        for n in self.innerNodes:
            if len(n.sequences) == len(self.sequences) - 1:
                l = len(n.pathLabel)
                if l >= minimumLength:
                    yield [(seq, idx, idx+l) for (seq,idx) in n.pathIndices]


def generate_overlap_graph(seqs, min_overlap=1):
    """
    Uses an all suffix-prefix search based on:

    Gusfield, D., Landau, G. M., & Schieber, B. (1992).
    An efficient algorithm for the all pairs suffix-prefix problem.
    Information Processing Letters, 41(4), 181-185.

    to construct an overlap graph of the input sequences

    :param seqs: List of string like objects
    :type: ~Bio.Seq
    :param int min_overlap: the minimum overlap needed to be connected in the overlap graph
    :return: Adjacency matrix of the overlap graph
    """
    def representStack(stacks):
        return "".join("["+",".join("Node({},{})".format(v.dfsRank,v.L) for v in s)+"]" for s in stacks)

    def updateBacktrace(n, past):
        # test if backward traversed and update the stacks accordingly

        if past.treeDepth > n.treeDepth != 0:
            n_parent = n.parent
            tmp = past.parent
            while n_parent != tmp:

                for j in tmp.L:
                    try:
                        stacks[j].pop()
                    except IndexError:
                        continue

                tmp = tmp.parent

    k = len(seqs)
    stacks = [[] for _ in xrange(k)]
    adja = np.zeros(shape=(k, k))
    st = GeneralisedSuffixTree(seqs)
    past = st.root
    reg = re.compile(r"(\D+)")

    for n in st.preOrderNodes:

        updateBacktrace(n, past)
        past = n

        if n.isLeaf and not n.suffix[0].isdigit():
            suffix = reg.search(n.suffix).group()
            if len(suffix) != len(seqs[n.sequences[0]]):
                continue

            # record suffixes of for prefix of string stored in SeqIdx
            prefix_idx = n.sequences[0]

            for suffix_idx, s in enumerate(stacks):
                if s:

                    m = s[0]
                    overlap = s[0].pathLabel

                    if len(overlap) >= min_overlap and suffix_idx != prefix_idx and prefix_idx < k:
                        adja[suffix_idx, prefix_idx] = len(seqs[prefix_idx]) - len(overlap)
                else:
                    adja[suffix_idx, prefix_idx] = len(seqs[prefix_idx])
        # forward traversal
        else:
            # push onto the ith stack, for each  L(n)
            for j in n.L:
                stacks[j].append(n)
            print "stacks:", representStack(stacks)

    return adja


import itertools
import numpy

INF = 1000000


class TSPExact:
    def __init__(self, matrix, redundant=False):
        """
        input: Adjacency matrix N x N
        len(nodes) = shape(input)
        @param matrix: N x N matrix
        @param redundant: false by default (n-1)!
        """
        # 1: Check if the adjacency matrix is symmetric
        self.cost_matrix = np.array(matrix)
        self.symmetry = None
        if self.symmetric() is False:
            self.symmetry = False
            self.as_symmetric()

        self.size = len(self.cost_matrix)
        self.nodes = range(self.size)
        self.all_tours = itertools.permutations

        if redundant is False:
            self.tours = self.non_redundant()
        else:
            self.tours = list(self.all_tours(self.nodes))

    def as_symmetric(self):
        """
        Reformulate an asymmetric TSP as a symmetric TSP:
        "Jonker and Volgenant 1983"
        This is possible by doubling the number of nodes. For each city a dummy
        node is added: (a, b, c) => (a, a', b, b', c, c')

        distance = "value"
        distance (for each pair of dummy nodes and pair of nodes is INF)
        distance (for each pair node and its dummy node is -INF)
        ------------------------------------------------------------------------
          |A'   |B'   |C'   |A    |B    |C    |
        A'|0    |INF  |INF  |-INF |dBA  |dCA  |
        B'|INF  |0    |INF  |dAB  |-INF |dCB  |
        C'|INF  |INF  |0    |dAC  |dBC  |-INF |
        A |-INF |dAB  |dAC  |0    |INF  |INF  |
        B |dBA  |-INF |dBC  |INF  |0    |INF  |
        C |dCA  |dCB  |-INF |INF  |INF  |0    |

        For large matrix an exact solution is infeasible
        if N > 5 (N = N*2 > 10) then use other techniques:
        Heuristics and relaxation methods
        @return: new symmetric matrix

        [INF][A.T]
        [A  ][INF]
        """

        shape = len(self.cost_matrix)
        mat = np.identity(shape) * - INF + self.cost_matrix

        new_shape = shape * 2
        new_matrix = np.ones((new_shape, new_shape)) * INF
        np.fill_diagonal(new_matrix, 0)

        # insert new matrices
        new_matrix[shape:new_shape, :shape] = mat
        new_matrix[:shape, shape:new_shape] = mat.T
        # new cost matrix after transformation
        self.cost_matrix = new_matrix

    def total_tour(self, tour):
        total = sum(self.cost_matrix[tour[node], tour[(node + 1) % self.size]] for node in self.nodes)
        return total, tour

    def shortest(self):
        min_tour = min(self.total_tour(tour) for tour in self.tours)
        if self.symmetry is False:
            min_tour_real = min_tour[0] + self.size / 2 * INF
            min_cycle_real = min_tour[1][0::2]
            return min_tour_real, min_cycle_real
        else:
            return min_tour

    def non_redundant(self):
        start_node = self.nodes[0]
        new_set = list(self.nodes)[1:]
        return [list(np.append(start_node, tour)) for tour in self.all_tours(new_set)]

    def symmetric(self):
        #@return: boolean symmetric | asymmetric
        mat = [tuple(row) for row in self.cost_matrix]
        return mat == zip(*mat)

    @property
    def tour_iterations(self):
        return list(self.tours)

    @property
    def tour_size(self):
        return self.size
