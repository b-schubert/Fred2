# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
   .. module:: Mosaic
   :synopsis:  This class implements the epitope selection functionality
    of to construct so called mosaic vaccines

    This module builds upon Coopr's Pyomo, an embedded algebraic modeling
    languages [2].

    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice


.. moduleauthor:: schubert

"""

from __future__ import division

import itertools as itr
import copy
import math

import numpy as np

import multiprocessing as mp

from Fred2.Core.Result import EpitopePredictionResult


def _popcover_score(args):
    """
    calculates the greedy PopCover score
    :param args: a tuple of arguments ()
    :return: the popcover score for peptide i
    """
    i, beta, p_h, p_g, E, R, n_allele, n_var = args

    return sum((R.get((h, g, i), 0)*p_g[g]*p_h[h])/(beta+E.get((h, g), 0))
               for h in xrange(n_allele)
               for g in xrange(n_var))


class PopCover(object):
    """
    Reimplementation of PopCover

    PopCover's immunogenicity function is defined as:


    S_j^{H,G} = \sum_{h \in H} \sum_{g \in G} \frac{R_{h,g}^j \codt p_h \codt p_g}{\beta + E_{h,g}}

    The sum is over all genomes g and HLA alleles h. Rjki is 1 if epitope j is present in genome g and is presented
    by allele h, and Ehg is the number of times allele h has been targeted by epitopes in genome g by the already
    selected set of epitopes, p_h is the frequency of allele h in a given population and p_g is the genomes frequency


    .. note::

    Lundegaard, C., Buggert, M., Karlsson, A. C., Lund, O., Perez, C., & Nielsen, M. (2010, August).
    PopCover: a method for selecting of peptides with optimal population and pathogen coverage.
    In Proceedings of the First ACM International Conference on Bioinformatics and Computational Biology
    (pp. 658-659). ACM.

    """

    def __init__(self, _results, threshold=None, beta=0.1, k=10, p_g=None, verbosity=0, processes=1):
        # check input data
        if not isinstance(_results, EpitopePredictionResult):
            raise ValueError("first input parameter is not of type EpitopePredictionResult")

        self.k = k

        # start constructing model
        _alleles = copy.deepcopy(_results.columns.values.tolist())

        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        if len(no_prob) > 0:
            #group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            for g, v in no_prob_grouped.iteritems():
                total_loc_a = len(v)
                if g in prob_grouped:
                    remaining_mass = 1.0 - sum(a.prob for a in prob_grouped[g])
                    for a in v:
                        a.prob = remaining_mass/total_loc_a
                else:
                    for a in v:
                        a.prob = 1.0/total_loc_a

        method = _results.index.values[0][1]
        res_df = _results.xs(_results.index.values[0][1], level="Method")
        res_df = res_df[res_df.apply(lambda x: any(x[a] > threshold.get(a.name, -float("inf"))
                                                   for a in res_df.columns), axis=1)]


        #construct structures needed for fast calculation
        variations = list(set(prot for p in res_df.index for prot in p.get_all_proteins()))
        v_to_idx = {p:i for i,p in enumerate(variations)}
        a_to_idx = {a:i for i,a in enumerate(_alleles)}
        if not variations:
            raise ValueError("Epitopes do not have protein origin annotated which is needed for algorithm.")

        peps = []
        alleles_epi = {}
        epi_alleles = {}
        epi_var = {}
        var_epi = {}
        R = mp.Manager().dict()
        for i,tup in enumerate(res_df.itertuples()):
            k = i
            p = tup[0]
            peps.append(p)
            for a, s in itr.izip(_alleles, tup[1:]):
                if method in ["smm", "smmpmbec", "arb", "comblibsidney"]:
                    try:
                        thr = min(1., max(0.0, 1.0 - math.log(threshold.get(a.name),
                                                              50000))) if a.name in threshold else -float("inf")
                    except:
                        thr = 0

                    if s >= thr:
                        a_idx = a_to_idx[a]
                        alleles_epi.setdefault(a_idx, set()).add(i)
                        epi_alleles.setdefault(i, set()).add(a_idx)
                        for pr in p.get_all_proteins():
                            pr_idx = v_to_idx[pr]
                            R[a_idx, pr_idx, k] = 1

                else:
                    if s >= threshold.get(a.name, -float("inf")):
                        a_idx = a_to_idx[a]
                        alleles_epi.setdefault(a_idx, set()).add(k)
                        epi_alleles.setdefault(k, set()).add(a_idx)
                        for pr in p.get_all_proteins():
                            pr_idx = v_to_idx[pr]
                            R[a_idx, pr_idx, k] = 1

            for pr in p.get_all_proteins():
                pr_idx = v_to_idx[pr]
                epi_var.setdefault(k, set()).add(pr_idx)
                var_epi.setdefault(pr_idx, set()).add(k)

        del v_to_idx
        del a_to_idx

        self.__verbosity = verbosity
        self.__processes = processes
        self.beta = beta

        self.__p_h = np.array([a.prob for a in _alleles])
        if p_g is None:
            self.__p_g = np.ones(shape=len(variations))
        else:
            raise NotImplementedError
        self.__peps = peps
        self.__alleles = _alleles
        self.__variaitons = variations
        self.__R = R
        self.__alleles_epi = alleles_epi
        self.__epi_alleles = epi_alleles
        self.__var_epi = var_epi
        self.__epi_var = epi_var

    def solve(self):
        """
        Uses a greedy approximation algorithm to construct a mosaic vaccine


        :return:
        """

        def __find_update_idx(i):
            """
            calculates the set of peptides who's score has to be updated

            :param int i: the index of the current selection
            :return: list of peptide indices for update
            """
            effected_alleles = self.__epi_alleles[i]
            effected_vars = self.__epi_var[i]

            epi_all = []
            for a in effected_alleles:
                epi_all.extend(self.__alleles_epi[a])

            epi_va = []
            for v in effected_vars:
                epi_va.extend(self.__var_epi[v])

            return list(set(epi_va).intersection(set(epi_all)) - set(selection_naive))

        pool = mp.Pool(processes=self.__processes)

        n_peps = len(self.__peps)
        n_alleles = len(self.__alleles)
        n_vars = len(self.__variaitons)
        R = self.__R
        selection_naive = []
        beta = self.beta
        scoring = np.zeros(shape=n_peps)
        p_h = self.__p_h
        p_g = self.__p_g
        E = mp.Manager().dict()

        # a list of peptides that need to be updated each round
        to_update = range(n_peps-1)

        #initialize the score
        scoring[to_update] = pool.map(_popcover_score,
                                      [(i, beta, p_h, p_g, E, R, n_alleles, n_vars) for i in
                                       to_update])
        curr_selection = np.argmax(scoring)
        selection_naive.append(curr_selection)
        scoring[curr_selection] = -float("inf")

        # simple selection
        while len(selection_naive) < self.k:
            # find indices that have to be updated
            for i in self.__epi_alleles[curr_selection]:
                for j in self.__epi_var[curr_selection]:
                    E[i, j] = E.get((i, j), 0) + 1

            to_update = __find_update_idx(curr_selection)
            scoring[to_update] = pool.map(_popcover_score,
                                          [(i, beta, p_h, p_g, E, R, n_alleles, n_vars) for i in
                                           to_update])
            curr_selection = np.argmax(scoring)
            scoring[curr_selection] = -float("inf")

            selection_naive.append(curr_selection)

        return [self.__peps[i] for i in selection_naive]