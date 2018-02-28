"""
Contains classes Quandle, Biquandle, Identity_Quandle,
Alexander_Quandle, and Singquandle.

FIXME:
    - Nothing for now.

TODO:
    - If X is a rack with operation a*b, then it is a birack if we
      define a**b as the identity a**b == a. Thus biquandle matrix2
      should be optional.
    - Does the above apply to singquandles/singracks?
    - homomorphism methods.
"""

import numpy as np
from pyknots.modules.magmas import Magma
from pyknots.modules.groups import Group
from pyknots.modules.utils import issquare, applymorphism
import json
import os

__all__ = ['Quandle', 'Biquandle', 'Singquandle',
'Alexander_Quandle', 'Conj_Quandle', 'Trivial_Quandle']

class Quandle(Magma):
    """Instantiate a quandle object by passing Quandle() a matrix, a
       string representation of the RIG index, or a numpy array.
    """
    def __init__(self, matrix):
        if type(matrix) is str:
            dirname = os.path.dirname(os.path.realpath(__file__)) 
            path = os.path.join(os.path.dirname(dirname), 'data', 'RIG_quandles.json')
            try:
                with open(path, 'r') as fp:
                    matrix = json.load(fp)[matrix]
                    self.__init__(matrix)
            except KeyError:
                raise TypeError('Input %s not a RIG quandle.' % (matrix))
        else:
            super().__init__(matrix)

    def __str__(self):
        return str(self.array)

    def as_biquandle(self):
        """ Generate identity matrix and return biquandle class object."""
        M1 = self.array
        M2 = Trivial_Quandle(self.order, self.index).array
        B = Biquandle(M1, M2)
        return B

    def inverse_quandle(self):
        M = self.array
        n = self.order
        new_M = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                k = M[i,j]
                new_M[k,j] = i
        return Quandle(new_M)

    def is_rack(self):
        """ A rack is a set with 2 axioms: for a, b, c in X,
            1) a*b is a bijection.
            2) (a*b)*c == (a*c)*(b*c).
        """
        M = self.array
        ind = self.index
        for i in range(self.order):
            col = []
            for j in range(self.order):
                for k in range(self.order):
                    if M[M[i,j]-ind,k] != M[M[i,k]-ind,M[j,k]-ind]:
                        return False
                col.append(M[j,i])
            for c in range(len(col)):
                if not c+ind in col:
                    return False
        return True

    def is_quandle(self):
        """ A quandle is a rack that satisfies the third axiom that for
            all a in X, a*a == a.
        """
        M = self.array
        ind = self.index
        if not self.is_rack():
            return False
        for i in range(self.order):
            if M[i,i] != i+ind:
                return False
        return True

    def is_biquandle(self):
        return False

    def is_singquandle(self):
        return False

    def is_kei(self):
        if self.is_quandle() and self.is_involutory():
            return True
        return False

    def is_trivial(self):
        """ A quandle is trivial if for all a, b in X, a*b == a."""
        M = self.array
        for i in range(self.order):
            for j in range(self.order):
                if M[i,j] != i:
                    return False
        return True

    def is_involutory(self):
        """ A quandle is involutory if for all a, b in X, a*(a*b) == b."""
        M = self.array
        for i in range(self.order):
            for j in range(self.order):
                if M[i,M[i,j]] != j:
                    return False
        return True

    def is_dihedral(self):
        """ If for all a, b in X, a*b == 2*b-a then it is dihedral.
            Equivalent to isomorphism to Alexander_Quandle(p, -1)
        """
        M = self.array
        ind = self.index
        p = self.order
        for i in range(self.order):
            for j in range(self.order):
                if M[i,j] != ((2*(j) - (i)) % p)+ind:
                    return False
        return True

    def is_medial(self):
        """ Equivalent to abelian quandle. If X satisfies the property that for
            any a, b, c, d in Q, (a*b)*(c*d) == (a*c)*(b*d) it is medial.
        """
        M = self.array
        ind = self.index
        for i in range(self.order):
            for j in range(self.order):
                for m in range(self.order):
                    for n in range(self.order):
                        if M[M[i,j]-ind,M[m,n]-ind] != M[M[i,m]-ind,M[j,n]-ind]:
                            return False
        return True


class Biquandle(object):
    """ Instantiate a biquandle object by passing Biquandle() a pair of
        matrices, a string representation of the RIG index, or a
        numpy array.
    """
    def __init__(self, matrix1, matrix2=None):
        if matrix2 is None:
            M1 = Quandle(matrix1)
            M2 = Identity_Quandle(M1.order, M1.index)
            self.__init__(M1.array, M2.array)
        else:
            self.quandle1, self.quandle2 = Quandle(matrix1), Quandle(matrix2)
            self.array1, self.array2 = np.array(matrix1), np.array(matrix2)
            self.order = len(matrix1[0])
            self.index = self._index()

    def __str__(self):
        return str(self.array1)+str(self.array2)

    def _index(self):
        """ Verify that indices of input match."""
        ind1, ind2 = np.amin(self.array1), np.amin(self.array2)
        if ind1 != ind2:
            raise IndexError('%s, %s have non-matching indices.' % (self.array1, self.array2))
        return ind1

    def is_birack(self):
        """ A birack is a set with 4 axioms: for a, b, c in X,
            1) a*b, a**b is a bijection.
            2) (a**b)**(c**b) == (a**c)**(b*c).
            3) (a*b)*(c*b) == (a*c)*(b**c)
            4) (a*b)**(c*b) == (a**c)*(b**c)
        """
        M1, M2 = self.array1, self.array2
        ind = self.index
        if not self.is_invertible():
            return False
        for a in range(self.order):
            for b in range(self.order):
                for c in range(self.order):
                    if M2[M2[a,b]-ind,M2[c,b]-ind] != M2[M2[a,c]-ind,M1[b,c]-ind]:
                        return False
                    if M1[M1[a,b]-ind,M1[c,b]-ind] != M1[M1[a,c]-ind,M2[b,c]-ind]:
                        return False
                    if M2[M1[a,b]-ind,M1[c,b]-ind] != M1[M2[a,c]-ind,M2[b,c]-ind]:
                        return False
        return True

    def is_biquandle(self):
        """ A biquandle is a birack such that for all a in X, there
            exists x such that: x*a == x <==> a**x == a"""
        M1, M2 = self.array1, self.array2
        if not self.is_birack():
            return False
        for i in range(self.order):
            for j in range(self.order):
                if M1[i,j] == i and M2[j,i] != j:
                    return False
        return True

    def is_singquandle(self):
        return False

    def is_invertible(self):
        if self.quandle1.is_left_invertible():
            if self.quandle2.is_left_invertible():
                return True
        return False

    def check_wada(self):
        M1, M2 = self.array1, self.array2

        for a in range(self.order):
            for b in range(self.order):
                for c in range(self.order):
                    if M1[M1[a,b],M1[M2[a,b],c]] != M1[a,M1[b,c]]:
                        return False
                    if M2[M1[a,b],M1[M2[a,b],c]] != M1[M2[a,M1[b,c]],M2[b,c]]:
                        return False
                    if M2[M2[a,b],c] != M2[M2[a,M1[b,c]],M2[b,c]]:
                        return False
        return True

class Singquandle(object):
    """ Instantiate a singquandle object by passing Singquandle() three
        matrices (denoting Cayley tables for *, R1, and R2), a string
        representation of the RIG index, or a numpy array.
        (Only matrices supported)
    """
    def __init__(self, matrix1, matrix2, matrix3):
        self.quandle1, self.quandle2 = Quandle(matrix1), Quandle(matrix2)
        self.quandle3 = Quandle(matrix3)
        self.array1, self.array2 = np.array(matrix1), np.array(matrix2)
        self.array3 = np.array(matrix3)
        self.order = len(matrix1[0])
        self.index = self._index()

    def __str__(self):
        return str(self.array1)+str(self.array2)+str(self.array3)

    def _index(self):
        """ Verify that indices of input match."""
        ind1, ind2, ind3 = np.amin(self.array1), np.amin(self.array2), np.amin(self.array3)
        if ind1 != ind2 != ind3:
            raise IndexError('%s, %s, %s have non-matching indices.' % (self.array1, self.array2, self.array3))
        return ind1

    def is_invertible(self):
        """ Check whether * is an invertible operation."""
        if not self.quandle1.is_left_invertible():
            return False
        return True

    """
    def check_identity(self):
         R1(x,y) = R2(y,x)*x, R2(x,y) = R1(y,x)*y.
        M1, M2, M3 = self.array1, self.array2, self.array3
        ind = self.index
        for a in range(self.order):
            for b in range(self.order):
                if M2[a,b] != M1[M3[b,a]-ind,a]:
                    return False
                if M3[a,b] != M1[M2[b,a]-ind,b]:
                    return False
        return True
    """

    def is_singquandle(self):
        """ Check if the object is a singquandle."""
        if self.is_nonoriented_singquandle() or self.is_oriented_singquandle():
            return True
        return False

    def is_nonoriented_singquandle(self):
        """ Check if the singquandle satisfies the axioms of a nonoriented
            singquandle.
        """
        M1, M2, M3 = self.array1, self.array2, self.array3
        ind = self.index
        if not self.is_invertible():
            return False
        for a in range(self.order):
            for b in range(self.order):
                if M3[a,b] != M2[b,M1[a,b]-ind]:
                    return False
                if M3[b,M1[a,b]-ind] != M1[M2[a,b]-ind,M3[a,b]-ind]:
                    return False
                if M2[a,b] != M3[M1[b,a]-ind,a]:
                    return False
                if M2[M1[b,a]-ind,a] != M1[M3[a,b]-ind,M2[a,b]-ind]:
                    return False
                for c in range(self.order):
                    if M1[M1[b,c]-ind,M3[a,c]-ind] != M1[M1[b,a]-ind,M2[a,c]-ind]:
                        return False
                    if M1[M2[a,b]-ind,c] != M2[M1[a,c]-ind,M1[b,c]-ind]:
                        return False
                    if M1[M3[a,b]-ind,c] != M3[M1[a,c]-ind,M1[b,c]-ind]:
                        return False
        return True

    def is_oriented_singquandle(self):
        """ Check if the singquandle satisfies the axioms of an oriented
            singquandle.
        """
        M1, M2, M3 = self.array1, self.array2, self.array3
        n, ind = self.order, self.index
        inv = self.quandle1.inverse_quandle().array
        if not self.is_invertible():
            return False
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if M1[M2[inv[x,y],z],y] != M2[x,M1[z,y]]:
                        return False
                    if M3[inv[x,y],z] != inv[M3[x,M1[z,y]], y]:
                        return False
                    if M1[inv[y,M2[x,z]],x] != inv[M1[y,M3[x,z]], z]:
                        return False
                    if M3[x,y] != M2[y,M1[x,y]]:
                        return False
                    if M1[M2[x,y], M3[x,y]] != M3[y, M1[x,y]]:
                        return False
        return True

class Alexander_Quandle(Quandle):
    """ Returns quandle generated by Alexander module Z_p/((t**b)-a).
        Setting exponent b not supported.
    """
    def __init__(self, p, a=-1, b=None):
        M = np.zeros((p, p), dtype=int)
        for i in range(p):
            for j in range(p):
                M[i,j] = (a*i + (1 - a)*j) % p
        super().__init__(M.tolist())

class Conj_Quandle(Quandle):
    """ Returns quandle generated by the group cayley table with automorphism f. Pass
        f as a permutation in the form (1, 2, 3, ...). Then f maps index(i) to i.
        Quandle is given by conjugation operation x*y = f(y)^{-1}f(x)f(y).
    """
    def __init__(self, matrix, *args):
        n = len(matrix[0])
        M = np.zeros((n, n), dtype=int)
        for arg in args:
            if isinstance(arg, tuple):
                matrix = applymorphism(matrix, arg)
        G = Group(matrix)
        m = G.array
        for i in range(n):
            for j in range(n):
                M[i,j] = m[m[G.inverse(j), i], j]
        super().__init__(M)


class Trivial_Quandle(Quandle):
    """ Returns a trivial quandle such that for all a, b in X,
        a*b == a. Optional index for non-0-indexed quandles.
    """
    def __init__(self, dim, index=0, flip=False):
        ind = index
        M = []
        if not flip:
            for i in range(dim):
                row = [i+ind for j in range(ind, dim+ind)]
                M.append(row)
        else:
            for i in range(dim):
                row = [j+ind for j in range(ind, dim+ind)]
                M.append(row)
        super().__init__(M)
