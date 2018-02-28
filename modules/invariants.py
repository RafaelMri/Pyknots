# Braid invariants.
from sympy import Poly
from sympy import primitive
from sympy import symbols
from sympy.abc import x, y, h, t
from sympy.matrices import Matrix, zeros
#from math import isnan
import numpy as np
#from itertools import product

__all__ = ['writhe', 'linking_number', 'seifert_matrix', 'alexander_matrix', 'alexander_polynomial',
            'alexander_system', ]

def writhe(braid):
    """ Find writhe of braid. Writhe = # of positive crossings minus
        # of negative crossings.
    """
    writhe = 0
    for i in braid.braid:
        if i < 0:
            writhe -= 1
        else:
            writhe += 1
    return writhe

def linking_number(braid):
    """ Linking number is given by writhe/2."""
    return writhe(braid)/2

def generators(braid):
    h = [0 for i in range(braid.length)]
    for i, j in enumerate(braid.braid):
        for k, m in enumerate(braid.braid[i+1:]):
            if abs(j) == abs(m):
                h[i] = i+k+1
                break
    return h

def seifert_matrix(braid):
    """ Computes Seifert matrix of the given braid. See [4] in the
        references for more details on the algorithm.
    """
    h = generators(braid)
    M = np.zeros((len(h), len(h)), dtype=int)
    b = braid.braid
    # check interaction between generators:
    for i in range(len(h)):
        if h[i] == 0:
            continue
        # chech each generator with itself:
        if b[i] > 0 and b[h[i]] > 0:
            M[i,i] = -1
        if b[i] < 0 and b[h[i]] < 0:
            M[i,i] = 1
        # check adjacent generators:
        if abs(b[i]) == abs(b[h[i]]) == abs(b[h[h[i]]]) and h[h[i]] != 0:
            if b[h[i]] < 0:
                M[i,h[i]] = -1
            else:
                M[h[i],i] = 1
        # check diagonal generators:
        for j in range(i+1, h[i]):
            if h[i] > h[j]:
                continue
            if abs(b[i]) - abs(b[j]) == 1:
                M[j,i] = -1
            if abs(b[i]) - abs(b[j]) == -1:
                M[i,j] = 1

    dels = [i for i in range(len(h)) if h[i] == 0]
    M = np.delete(M, dels, axis = 0)
    M = np.delete(M, dels, axis = 1)
    return M

def _normalize_apoly(apoly):
    """ Normalizes alexander polynomial. Factors out greatest
        power of t possible and multiplies by -1 if necessary.
    """
    apoly = primitive(Poly(apoly, t))[-1]
    apoly = apoly.terms_gcd()[-1]
    if apoly.EC() < 0:
        apoly = -1 * apoly
    return apoly.as_expr()

def alexander_matrix(braid, t=t, modulus=None):
    """ Returns presentation matrix of alexander polynomial."""
    # Specifying t and modulus not working. See class docstring for info
    sys = alexander_system(braid, t, modulus)
    all_vars = ['z'+str(i) for i in range(len(braid))]
    M = zeros(len(braid), len(braid))

    for i in range(len(braid)):
        for j in range(len(braid)):
            if modulus is None:
                M[i,j] = sys[i].coeff(all_vars[j])
            else:
                M[i,j] = sys[i].coeff(all_vars[j]) % modulus

    return M

def alexander_polynomial(braid, t=t, modulus=None, seifert=False):
    """ Calculate Alexander Polynomial."""
    if seifert:
        M = Matrix(seifert_matrix(braid).tolist())
        if modulus is None:
            apoly = M - t*M.transpose()
        else:
            apoly = (M - t*M.transpose()) % modulus
        apoly = apoly.det()
    else:
        deleted = False
        h = generators(braid)
        M = alexander_matrix(braid, t, modulus)
        for i in range(len(braid)):
            if h[i] == 0:
                deleted = True
                M.row_del(i)
                M.col_del(i)
                break
        if not deleted:
            M.row_del(-1)
            M.col_del(-1)
        apoly = M.det()
    return _normalize_apoly(apoly)

def alexander_system(braid, t=t, modulus=None):
    """ Generates a system of equations for computing the alexander
        polynomial.
    """
    start_vars = braid.start_vars()
    strands = start_vars[:]
    sys, reduced_sys = [], []
    #all_vars = start_vars[:]

    def apoly(x,y,h,t):
        if modulus is None:
            return t*x + (1-t)*y - h
        else:
            return t*(x % modulus) + (1-t)*(y % modulus) - h

    for i, j in enumerate(braid.braid):
        z = 'z'+str(len(strands) + i+1)
        z = symbols(z)
        s1, s2 = strands[j-1], strands[j]
        if j > 0:
            strands[j-1], strands[j] = s2, z
            sys.append(apoly(s1, s2, z, t))
        elif j < 0:
            strands[j-1], strands[j] = z, s1
            sys.append(apoly(z, s1, s2, t))
    for i in sys:
        reduced_sys.append(i.subs(zip(strands, start_vars)))
    return reduced_sys

def determinant(braid):
    return alexander_polynomial(braid, -1)
