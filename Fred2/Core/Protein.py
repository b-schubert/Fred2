# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Protein
   :synopsis: Contains the Protein Class.
.. moduleauthor:: brachvogel

"""
__author__ = 'brachvogel'

import itertools as iter
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger




class Protein(MetadataLogger, Seq):
    """
    Protein corresponding to exactly one transcript.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param str gene_id: ID of the genome
    :param str transcript_id: ID of the corresponding transcript 
    :param Transcript orig_transcript: Reference to the originating transcript
    :param dict(int,list(Variant)) _vars: Nonsynonymous variants that are
                                          assoziated with the protein. 
                                          key=position within protein, 
                                          value=list of variants at that pos
    """

    def __init__(self, _seq, _gene_id, _transcript_id, _orig_transcript=None, 
                 _vars=None):
        """
        :param str _seq: String of an IUPACProtein alphabet, representing the
                         protein
        :param str _gene_id: ID of the genome the protein originated from
        :param str _transcript_id: ID of the transcript the protein originated 
                                   from
        :param Transcript _orig_transcript: Reference to the originating 
                                            transcript
        :param dict(int,list(Variant)) _vars: Nonsynonymous variants that are
                                              assoziated with the protein. 
                                              key=position within protein, 
                                              value=list of variants at that pos
        """
        # Init parent type:
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        # Init own member:
        if _vars is None:
            self.vars = dict()
        else:
            self.vars = _vars  # {prot-position: list(variant)}
        self.orig_transcript = _orig_transcript
        self.transcript_id = _transcript_id
        self.gene_id = _gene_id


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        :param int index: position within the primary sequence
        :returns: Protein -- A protein consisting of the single letter at
                  position :attr:`index`.
        """
        def frameshift_influences(res, start):
            # find variants out side the peptide frame, still influencing it via a
            # frameshift
            accu = {} # accumulator for relevant variants

            vars = sorted(iter.chain(*self.vars.values()), key=lambda v: v.genomePos)
            shift = 0

            for var in vars:

                pos = var.get_protein_position(self.transcript_id)
                new_shift = var.get_shift()

                if pos < start:
                    # does a variant yield a frame shift?
                    if shift + new_shift:
                        shift += new_shift
                        accu.setdefault(pos, []).append(var)
                    else:
                        accu = {}
                # here: var.get_protein_position >= start, we are done!
                else:
                    res.update(accu)
                    break

        if not isinstance(index, slice):
            index = slice(index, index+1)
        item = str(self)[index]
        tmp_var = {pos:var_list for pos, var_list in self.vars.iteritems() if index.start <= pos <= index.stop}

        frameshift_influences(tmp_var, index.start)
        p= Protein(item, self.gene_id, self.transcript_id, _orig_transcript=self.orig_transcript, _vars=tmp_var)
        return p
    def __repr__(self):
        # Header:
        lines = []
        lines += ["PROTEIN: %s (aa-seq)" % str(self)]
        lines += ["\t  %s (orig transcript)"%self.transcript_id]

        # Variants:

        lines += ["\t VARIANTS:"]
        for vpos, vset in self.vars.iteritems():
            for v in vset:
                lines.append('\t pos %i: %s'%(vpos, v))
        return '\n\t'.join(lines) + '\n'

    def __eq__(self, other):
        return str(self) == str(other)

    def __cmp__(self, other):
        return cmp(str(self), str(other))

    def __hash__(self):
        return hash(str(self))

