import os
import subprocess
from tempfile import NamedTemporaryFile

def solve_TSP_LKH(A, startNode=None, warmstart=None):
    """
    Solving a TSP problem by applying Lin-Kernighan traveling salesman heuristic

    .. note::

        LKH implementation must be downloaded, compiled, and globally executable.
        Source code can be found here:
        http://www.akira.ruc.dk/~keld/research/LKH/
    :param numpy.matrix A: A adjacency matrix
    :param int startNode: index of start node if fixed
    :param List(int) warmstart: a list of indixes defining a tour
    :return: A list of indixes representing the order of the hamiltonian circle
    :rtype: list(int)
    """

    tmp_conf = NamedTemporaryFile(delete=False)
    tmp_prob = NamedTemporaryFile(delete=False)
    tmp_out = NamedTemporaryFile(delete=False)
    n,m = A.shape

    # write config file:
    if warmstart is not None:
        tmp_warmstart = NamedTemporaryFile(delete=False)
        warmstart_file = "INITIAL_TOUR_FILE = %s\n"%tmp_warmstart.name
        with tmp_warmstart as f:
            f.write("\n".join(map(lambda x: str(x),warmstart))+"\n-1")
    else:
        warmstart_file = ""
    tmp_conf.write("PROBLEM_FILE = %s\n%sOUTPUT_TOUR_FILE = %s\n" %(tmp_prob.name, warmstart_file, tmp_out.name))
    tmp_conf.close()
    epis = []
    # write problem file:
    tmp_prob.write \
        ( "NAME: %s\nTYPE: ATSP\nDIMENSION: %i\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION\n" %
        (tmp_prob.name, n))
    for i in xrange(n):
        epis.append(i)
        tmp_prob.write("\t".join("99999999" if i == j else str(int(A[i, j]*1000)) for j in xrange(n)) +"\n")

    tmp_prob.write("EOF")
    tmp_prob.close()

    # try:
    r = subprocess.call( "LKH %s" %tmp_conf.name, shell=True)
    if r == 127:
        raise RuntimeError("LKH is not installed or globally executable.")
    elif r != 0:
        raise RuntimeError("An unknown error occurred for method LKH. "
                           "Please check whether LKH is globally executable.")
    result = []
    with open(tmp_out.name, "r") as resul_f:
        is_tour = False
        for l in resul_f:
            if is_tour:
                if int(l.strip()) == -1:
                    break
                result.append(int(l.strip())-1)
            elif "TOUR_SECTION" in l:
                is_tour = True
            else:
                pass

    tmp_out.close()
    os.remove(tmp_conf.name)
    os.remove(tmp_prob.name)
    os.remove(tmp_out.name)
    return result