# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Generator
   :synopsis: Contains functions to transform variants to transcripts and 
              proteins to peptides. Transition of transcripts to proteins
              is done via :meth:`~Fred2.Core.Transcript.translate`
.. moduleauthor:: schubert

"""
__author__ = 'schubert'

import warnings
import math
import itertools as iter
from collections import defaultdict

from Fred2.Core.Base import COMPLEMENT
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Transcript import Transcript
from Fred2.Core.Variant import VariationType
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields

################################################################################
#Private module functions. It should not be possible to import these!

def _update_var_offset(vars, transId_old, transId_new):
    """
    :param var:
    :param transId_old:
    :param transId_new:
    """
    for v in vars:
        offset = v.offsets[transId_old]
        v.offsets[transId_new] = offset

def _incorp_snp(seq, var, transId, offset, externalOff=0):
    """
    incorporates a snp into the given transcript sequence

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: (list,int) the modified seq, the modified offset
    """
    #tmp_seq = seq[:]
    if VariationType.SNP != var.type:
        raise TypeError("%s is not a SNP"%str(var))
    var.offsets[transId] = offset

    #print transId, len(seq), var.get_transcript_position(transId)-1
    #if seq[var.get_transcript_position(transId)-1-externalOff] != var.ref:
    #    warnings.warn("For %s bp dos not mmatch ref of assigned variant %s. Pos %i, var ref %s, seq ref %s " % (
    #    transId, str(var), var.get_transcript_position(transId) - 1-externalOff, var.ref,
    #    seq[var.get_transcript_position(transId) - 1-externalOff]))

    #print "incopr SNP ", var,var.get_transcript_position(transId),var.get_transcript_position(transId)-1-externalOff,externalOff
    seq[var.get_transcript_position(transId)-1-externalOff] = var.obs

    return seq, offset


def _incorp_insertion(seq, var, transId, offset, externalOff=0):
    """
    incorporates an insertion into the given transcript sequence

    :param list(char) seq: (list) transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: (list,int) modified sequence, the modified offset
    """
    #tmp_seq = seq[:]
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    var.offsets[transId] = offset
    pos = var.get_transcript_position(transId)-externalOff

    seq[pos:pos] = var.obs
    return seq, offset + len(var.obs)


def _incorp_deletion(seq, var, transId, offset, externalOff=0):
    """
    incorporates a deletion into the given transcript sequence

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript 
                       position of the variant
    :return: (list,int) -- modified sequence, the modified offset
    """
    #tmp_seq=seq[:]
    if var.type not in [VariationType.DEL, VariationType.FSDEL]:
        raise TypeError("%s is not a deletion"%str(var))

    var.offsets[transId] = offset
    pos = var.get_transcript_position(transId)
    s = slice(pos-externalOff, pos+len(var.ref)-externalOff)
    del seq[s]
    return seq, offset - len(var.ref)


_incorp = {
            VariationType.DEL: _incorp_deletion,
            VariationType.FSDEL: _incorp_deletion,
            VariationType.INS: _incorp_insertion,
            VariationType.FSINS: _incorp_insertion,
            VariationType.SNP: _incorp_snp
         }


def _check_for_problematic_variants(vars):
    """
    Filters problematic variants, e.g. variants that coincide.

    :param list(Variant) vars: initial list of variants
    :return: boole -- ture if now intersecting variants were found
    :invariant: list(Variant) vars: List is sorted based on genome position
    """
    v_tmp = vars[:]
    v = v_tmp.pop()
    current_range = (v.genomePos, v.genomePos
                                      +len(v.ref)-1 if v.type in [VariationType.FSDEL, VariationType.DEL] else
                                      v.genomePos)
    for v in reversed(v_tmp):
        if v.genomePos <= current_range[1]:
            print "crash",current_range, v
            return False
        else:
            current_range = (v.genomePos, v.genomePos
                                      +len(v.ref)-1 if v.type in [VariationType.FSDEL, VariationType.DEL] else
                                      v.genomePos)
            #print "new block",v, current_range
    return True


#################################################################
# Public transcript generator functions

################################################################################
#        V A R I A N T S     = = >    P E P T I D E S
################################################################################
def generate_peptides_from_variants(vars, length, dbadapter, trans_filter=None):
    """
    generates polymorphic peptides based on given variants

    :param list(Variant) vars: list of Variants
    :param int length: length of peptides
    :param dbadapter: DBAdapter to fetch the transcript sequences
    :return: List(Peptide) -- A list of polymorphic peptides
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_peptides_from_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            v = vs.pop()
            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                # update the transcript variant id
                old_trans = generate_peptides_from_variants.transOff
                generate_peptides_from_variants.transOff += 1
                transOff = generate_peptides_from_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    def _generate_heterozygous(tId, vs, seq, usedVs, newUsedVs, externalOff, oldTiD, offset=0):
        """
            incorporates heterozygous variants into polymorphic sequences
        """
        #print "In generate_herto",tId
        if vs:
            v = vs.pop()
            #find offset
            if not offset:
                for uv in usedVs:
                    if uv.genomePos < v.genomePos:
                        offset = uv.offsets.get(tId, 0)
                    else:
                        break

            #generate combinatorial branches:
            vs_tmp = vs[:]
            tmp_seq = seq[:]
            tmp_newUsedVs = newUsedVs[:]

            for s in _generate_heterozygous(tId, vs_tmp, tmp_seq, usedVs,tmp_newUsedVs, externalOff, oldTiD, offset=offset):
                    yield s

            #old_trans = generate_peptides_from_variants.transOff

            generate_peptides_from_variants.transOff_hetero += 1
            transOff = generate_peptides_from_variants.transOff_hetero
            new_tId = "_".join(tId.split("_")[:-1])+"_%i"%transOff
            _update_var_offset(usedVs,  oldTiD, new_tId)

            seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v,  "_".join(tId.split("_")[:-1])+"_%i"%transOff, offset,
                                                                         externalOff=externalOff)

            newUsedVs.append(v)
            for s in _generate_heterozygous(tId, vs, seq, usedVs,newUsedVs, externalOff,  new_tId, offset=offset):
                yield s
        else:
            yield oldTiD, seq, newUsedVs
    ####################################################################################################################

    def frameshift_influences(tid, _vars, res, start):
        # find variants out side the peptide frame, still influencing it via a
        # frameshift
        accu = [] # accumulator for relevant variants

        _vars.sort(key=lambda v: v.genomePos) # necessary?
        shift = 0

        for var in _vars:

            pos = var.get_protein_position(tid)
            new_shift = var.get_shift()

            if pos < start:
                # does a variant yield a frame shift?
                if shift + new_shift:
                    shift += new_shift
                    accu.append(var)
                else:
                    accu = {}
            # here: var.get_protein_position >= start, we are done!
            else:
                res += accu
                break

    def get_influencing_vars(tId, vars, i, end):
        # Generate peptide sequences and find the variants within each

        pep_var = [v for v in vars if i <= v.get_transcript_position(tId) <= end]

        # outside variants that affect the peptide via frameshift:
        frameshift_influences(tId, vars, pep_var, i)
        return pep_var

    def trans_to_prot_pos(i):
        return int(max(math.ceil(i/3.0),1))

    ####################################################################################################################

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        #print v, v.isHomozygous
        for trans_id in v.coding.iterkeys():
            #print trans_id
            if len(trans_filter):
                if trans_id in trans_filter:
                    transToVar.setdefault(trans_id, set()).add(v)
            else:
                transToVar.setdefault(trans_id, set()).add(v)


    print "Nof transcripts",len(transToVar)
    #print transToVar
    prots = []
    c = 1
    for tId, vs in transToVar.iteritems():
       # print tId, c, len(transToVar)
        #print "VS total", len(vs)

        #print "Nof Transcripts ",len(transToVar)
        #print vs_hetero
        c += 1
        vs = list(vs)

        vs.sort(key=lambda v: v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue

        query = dbadapter.get_transcript_information(tId)
        #print "end fetching"
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for v in vs:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        generate_peptides_from_variants.transOff = 0

        vs_homo_and_fs = filter(lambda x: x.type in [VariationType.FSINS, VariationType.FSDEL] or x.isHomozygous, vs)
        generate_peptides_from_variants.transOff_hetero = 2**sum(1 for x in vs if x.type in [VariationType.FSINS, VariationType.FSDEL] and not x.isHomozygous)

        vs_hetero = filter(lambda x: not x.isHomozygous and x.type not in [VariationType.FSINS, VariationType.FSDEL], vs)
        #print "homo",len(vs_homo_and_fs)
        #print "hetero snp",len(vs_hetero)
        #print "VS total", len(vs)
        #print "homo",len(vs_homo_and_fs)
        #print "hetero snp",len(vs_hetero)
        #print "Nof Transcripts ",len(transToVar)
        #print vs_hetero
        for tId, varSeq, varComb in _generate_combinations(tId, vs_homo_and_fs, list(tSeq), [], 0):
            #print "In Loop",tId,varComb
            if vs_hetero:
                for i in xrange(0,len(varSeq)-(3*length-3),3):
                    end = i+3*length
                    tmp_var = varSeq[:]
                    #print "In Peptide loop",i,end, len(tmp_var),len(varSeq)+3*length-(3*length)
                    frac_var = filter(lambda x: i <= x.get_transcript_position(tId)-1 <= end, vs_hetero)

                    frac_usedVar = get_influencing_vars(tId, varComb, i, end+1)
                    if not frac_var and not frac_usedVar:
                        continue
                    #print "Fraction of heterozygious in rage ",i," ", end," ", [x.get_transcript_position(tId) for x in frac_var]

                    frac_usedVar = get_influencing_vars(tId, varComb, i, end)
                    for ttId, vvarSeq, vvarComb in _generate_heterozygous(tId, frac_var, tmp_var, frac_usedVar, [], 0,tId[:]):
                        #print "In Peptide combinatoric loop",ttId, frac_usedVar, vvarComb,"".join(vvarSeq),tId
                        vvarComb.extend(frac_usedVar)
                        #print vvarComb
                        vvarComb.sort(key=lambda x: x.genomePos)

                        if vvarComb:
                            try:
                                p = Transcript(geneid, ttId, "".join(vvarSeq), _vars=vvarComb).translate()
                                #print "Protein ",repr(p),i,end,i//3,i//3+length,p[i//3:i//3+length]
                                prots.append(p[i//3:i//3+length])
                            except ValueError:
                                #print "Protein not multiple of 3", ttId, vvarComb
                                continue
            else:
                try:
                    prots.append(Transcript(geneid, tId, "".join(varSeq), _vars=varComb).translate())
                except ValueError:
                    continue
        #print tId, len(prots)

    #return generate_peptides_from_protein(prots, length)
    return filter(lambda x: "*" not in str(x), generate_peptides_from_protein(prots, length))

################################################################################
#        V A R I A N T S     = = >    T R A N S C R I P T S
################################################################################


def generate_transcripts_from_variants(vars, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param list(Variant) vars: A list of variants for which transcripts should 
                               be build
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: (Generator(Transcripts)) -- a generator of transcripts with all 
             possible variations determined by the given
             variant list
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_transcripts_from_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            v = vs.pop()
            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                # update the transcript variant id
                old_trans = generate_transcripts_from_variants.transOff
                generate_transcripts_from_variants.transOff += 1
                transOff = generate_transcripts_from_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        #print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_transcripts_from_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), [], 0):
            yield Transcript(geneid, tId, "".join(varSeq), _vars=varComb)


def generate_transcripts_from_tumor_variants(normal, tumor, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param list(Variant) normal: A list of variants of the normal tissue
    :param list(Variant) tumor: A list of variant of the cancer tissue for which transcript should be generated
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: (Generator(Transcripts)) -- a generator of transcripts with all
             possible variations determined by the given
             variant list
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_transcripts_from_tumor_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            flag, v = vs.pop()

            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                if flag:
                    for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                        yield s

                # update the transcript variant id
                old_trans = generate_transcripts_from_tumor_variants.transOff
                generate_transcripts_from_tumor_variants.transOff += 1
                transOff = generate_transcripts_from_tumor_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in tumor:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append((False, v))

    for v in normal:
        for trans_id in v.coding.iterkeys():
            if trans_id in transToVar:
                transToVar.setdefault(trans_id, []).append((True, v))

    for tId, vs in transToVar.iteritems():
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for flag, v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v[1].genomePos, reverse=True)
        if not _check_for_problematic_variants(map(lambda x: x[1],vs)):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue

        generate_transcripts_from_tumor_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), [], 0):
            yield Transcript(geneid, tId, "".join(varSeq), _vars=varComb)


################################################################################
#        P R O T E I N    = = >    P E P T I D E
################################################################################

def generate_peptides_from_protein(proteins, window_size):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.

    :param Protein protein: list of proteins from which a list of unique
                            peptides should be generated
    :param int window_size: size of peptide fragments
    """
    def isValidAASequence(epitope):

        epitope = epitope.strip()
        if epitope == '':
            return 1

        if not epitope.isalpha():
            return 0

        # check for invalid letters
        for aa in ['B','J','O','U','X','Z']:
            if aa in epitope:
                return 0
        return 1

    def frameshift_influences(tid, _vars, res, start):
        # find variants out side the peptide frame, still influencing it via a
        # frameshift
        accu = [] # accumulator for relevant variants

        _vars.sort(key=lambda v: v.genomePos) # necessary?
        shift = 0

        for var in _vars:

            pos = var.get_protein_position(tid)
            new_shift = var.get_shift()

            if pos < start:
                # does a variant yield a frame shift?
                if shift + new_shift:
                    shift += new_shift
                    accu.append(var)
                else:
                    accu = {}
            # here: var.get_protein_position >= start, we are done!
            else:
                res += accu
                break


    def gen_peptide_info(protein):
        # Generate peptide sequences and find the variants within each
        res = []

        seq = str(protein)
        for i in xrange(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]

             # get the variants affecting the peptide:
            if protein.vars:
                # variants within the peptide:
                pep_var = [var for pos, var_list in protein.vars.iteritems() \
                           for var in var_list if i <= pos <= end+1]
                # outside variants that affect the peptide via frameshift:
                frameshift_influences(protein.transcript_id, 
                                      protein.orig_transcript.vars.values(),
                                      pep_var, i)
            else:
                pep_var = []

            res.append((pep_seq, pep_var))
        return res


    final_peptides = {} # sequence : peptide-instance

    for prot in proteins:
        # generate all peptide sequences per protein:
        for (seq, _vars) in gen_peptide_info(prot):
            if not isValidAASequence(seq):
                continue
                
            t_id = prot.transcript_id
            if seq not in final_peptides:
                final_peptides[seq] = Peptide(seq)

            final_peptides[seq].proteins[t_id] = prot
            final_peptides[seq].vars[t_id] = _vars
            final_peptides[seq].transcripts[t_id] = prot.orig_transcript

    return final_peptides.values()
