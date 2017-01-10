from pyknots import *
import numpy as np
from numpy.linalg import inv
from itertools import product

__all__ = ['singify', 'genquandles', 'permtomat', 'issquare', 'allelems', 'mergedicts',
           'applymorphism', 'directproduct']

def _check_type(*args):
    """ Checks type of input, returns tuple of input as np.ndarrays. Used
        to throw exceptions on incorrect input.
    """
    Mats = []
    for M in args:
        if isinstance(M, list):
            Mats.append(np.array(M))
        elif isinstance(M, np.ndarray):
            Mats.append(M)
        else:
            raise TypeError('Input %s not a list, np.ndarray, or Quandle object.' % (type(M)))

def singify(matrix1, matrix2, qtype=0):
    """ Given Cayley table for 2 ops, return the third.
        type:
          0  =>  input=(*, R1),     output=R2
          1  =>  input=(*, R2),     output=R1
          input=(R1, R2) and output=* not well defined.
    """
    _check_type(matrix1, matrix2)
    M1, M2 = np.array(matrix1), np.array(matrix2)
    order = len(M1[0])
    M3 = np.zeros((order, order), dtype=int)
    ind = abs(np.amax(M1) - (order - 1))
    if qtype == 0:
        for a in range(order):
            for b in range(order):
                M3[a,b] = M2[b,M1[a,b]-ind]
    elif qtype == 1:
        for a in range(order):
            for b in range(order):
                M3[a,b] = M2[M1[b,a]-ind,a]
    return M3

def genquandles(order, **kwargs):
    """ Return generator of all quandles of specified order. Possible kwargs
        include dihedral, connected, etc."""
    all_quandles = []
    perms = []
    M = np.zeros((order, order), dtype=int)
    """
    for i in range(order):
        M[i,i] = i
    for i in range(order):
        for j in range(order):
            for k in range(order):
                M[j,k] = i
                M[i,j]
                """
    raise NotImplemented

def permtomat(sigma):
    """ Constructs zero matrix representation of the
        permutations. Allows simpler permutation composition.
    """
    p = np.zeros((len(sigma), len(sigma)), dtype=int)
    for i in range(len(sigma)):
        p[i,sigma[i]] = 1
    return p

def issquare(*args):
    """ Check whether matrices are square. Input either list or np.ndarray."""
    for M in args:
        _check_type(M)
        rows = len(M)
        for row in M:
            if len(row) != rows:
                return False
    return True

def allelems(matrix):
    """ Returns all elements of matrix. Accepts input as list or np.ndarray."""
    #_check_type(matrix)
    elems = set([i for j in matrix for i in j])
    return list(elems)

def mergedicts(dicts):
    """ Merges dictionarys used for graph representations. Matching keys
        combines the lists. i.e: 1:[0] + 1:[1, 2] = 1:[0,1,2]. Mulitiplicity
        not allowed. Input is list of dicts.
    """
    new = {i: [] for i in dicts[0].keys()}
    for d in dicts:
        for i in d.keys():
            new[i].append(d[i])
    new = {i: allelems(new[i]) for i in new.keys()}
    return new

def applymorphism(matrix, perm):
    """ Generate isomorphic tables using specified permutation.
        Swaps matrix elements according to perm and returns:
        p**(-1) * M1 * p == M2. Input perm as a tuple, list, or np.ndarray, ie:
        p = (0, 1, 2) == [[1,0,0],[0,1,0],[0,0,1]] would be the
        identity permutation.
    """
    if not isinstance(perm, tuple):
        p = perm.copy()
        perm = tuple([list(i).index(1) for i in perm])
    else:
        p = permtomat(perm)
    ind = np.amin(matrix)
    vfunc = np.vectorize(lambda x: perm[x-ind]+ind)
    M = vfunc(matrix)
    M = np.dot(np.dot(inv(p), M), p)
    return M.astype(int)

def directproduct(group1, group2):
    """ Returns table representing the direct product of group1 and group2."""
    M1, M2 = group1.array, group2.array
    p = list(product(allelems(group1.array), allelems(group2.array)))
    order = group1.order * group2.order
    M = np.zeros((order, order), dtype=int)
    for i, j in enumerate(p):
        for m, n in enumerate(p):
            a = (M1[j[0], n[0]], M2[j[1], n[1]])
            M[i,m] = p.index(a)
    return M

def graphtomat(graph):
    """ Convert dictionary representations of graphs like those used in graph
        module to matrices, using 0 to fill in empty or unknown elements.
    """
    order = len(graph.keys())
    M = np.zeros((order, order), dtype=int)
    for i in graph.keys():
        for j in graph[i]:
            raise NotImplemented
