""" Group operation table generator

FIXME:
    - Define separate subgroup() for Dihedral_Group and Symmetric_Group
      that preserves perms attribute


TODO:
    - define Group.__mul__() for behavior:
      Z_3 * Z_3, D_3 * S_9 etc.
    - define Group.__div__() likewise, for quotients.
    - G.cosets(), G.all_subgroups(), G.is_subgroup(S), etc.
    - Add arg/kwarg support all over, ie:
      inverse(x, left=False)

"""

import numpy as np
from itertools import permutations
from math import factorial
from pyknots.modules.magmas import Magma
from pyknots.modules.utils import permtomat, allelems


__all__ = ['Group', 'Cyclic_Group', 'Multiplicative_Group',
'Dihedral_Group', 'Symmetric_Group']

class Group(Magma):
    """ Pass finite set as tuple, set, or list of elements, along with
        op - any function that accepts x, y, and modulus as args and
        returns 1 element. Will accept a list of lists and no op to
        directly pass a matrix.
    """
    def __init__(self, finite_set, op=None):
        if op is None:
            super().__init__(finite_set)
            self.finite_set = allelems(finite_set)
            self.identity = self._identity()
            self.op = op
        elif callable(op):
            self.finite_set, self.op = list(finite_set), op
            self.order = len(self.finite_set)
            self.index = min(self.finite_set)
            super().__init__(self._matrix())
            self.identity = self._identity()
        else:
            raise TypeError('Input %s not a function.' % (op))

    def _matrix(self):
        """ Set matrix representation of the group operation table."""
        M = np.zeros((self.order, self.order), dtype=int)
        p = self.order + self.index
        for i, j in enumerate(self.finite_set):
            for m, n in enumerate(self.finite_set):
                M[i,m] = self.op(j,n,p)
        return M

    def _identity(self):
        """ Finds and returns the identity element of the group."""
        M = self.array
        ind = self.index
        for i in range(self.order):
            identity = True
            for j in range(self.order):
                if M[i,j] != M[j,i]:
                    identity = False
                    break
                if M[i,j] != j+ind:
                    identity = False
                    break
            if identity == True:
                return i+ind
        return None

    def is_associative(self):
        """ For is_group() method, checks for associativity."""
        M = self.array
        ind = self.index
        for i in range(self.order):
            for j in range(self.order):
                for k in range(self.order):
                    if M[M[i,j]-ind,k] != M[i,M[j,k]-ind]:
                        return False
        return True

    def has_identity(self):
        """ For is_group() method, checks existence of identity element."""
        if self.identity == None:
            return False
        return True

    def has_inverses(self):
        """ For is_group() method, finds inverse for all elements."""
        M = self.array
        ind = self.index
        if self.identity == None:
            return False
        for i in range(self.order):
            inverse = False
            if i != self.identity - ind:
                for j in range(self.order):
                    if M[i,j] == self.identity:
                        inverse = True
                        break
                if inverse == False:
                    return False
        return True

    def is_group(self):
        """ Checks the group axioms: associativity on the set,
            existence of an identity element, and inverses for all
            elements.
        """
        if self.is_associative():
            if self.has_identity():
                if self.has_inverses():
                    return True
        return False

    def is_subgroup(self, other, lazy=True):
        """ Checks if other is subgroup of self. If lazy=False,
            the subgroup is normalized and the group axioms are
            checked. Cannot differentiate between group operations.
        """
        ind1, ind2 = self.index, other.index
        M1, M2 = other.array, self.array
        gens = other.finite_set
        norm = []
        for i in gens:
            if i not in self.finite_set:
                return False
        if len(gens) != other.order:
            return False
        for i, j in enumerate(gens):
            for m, n in enumerate(gens):
                print(i,m, j,n)
                if M1[i,m] != M2[j,n]:
                    return False
        if not lazy:
            # Normalize subgroup before checking group axioms.
            # Ie: [[0,4],[4,0]] => [[0,1],[1,0]]
            for i in range(self.order):
                n = [self.finite_set.index(j) for j in gens]
                norm.append(n)
            return Group(norm).is_group()
        return True

    def subgroup(self, gens, lazy=True):
        """ Returns subgroup generated by gens or None.
            If lazy is False then it finds a subgroup containing
            the gens as well as any other necessary elements.
        """
        M1 = np.zeros((len(gens), len(gens)), dtype=int)
        M2 = self.array
        elems = []
        for i in gens:
            if i not in self.finite_set:
                return None
        for i, j in enumerate(gens):
            for m, n in enumerate(gens):
                k = M2[j,n]
                M1[i,m] = k
                if k not in elems:
                    elems.append(k)
        if len(gens) == len(elems):
            return M1
        elif not lazy:
            return self.subgroup(elems, lazy=False)
        return None

    def is_abelian(self):
        """ If the operation is commutative on the finite set,
            then the group is abelian.
        """
        M = self.array
        for i in range(self.order):
            for j in range(self.order):
                if M[i,j] != M[j,i]:
                    return False
        return True

    def inverse(self, x):
        """ Find left or right inverse of x in the group. The returned
            inverse is not necessarily both left and right inverse.
        """
        # Might be better to allow the user to set left or
        # right inverses strictly with kwargs.
        ind = self.index
        M = self.array
        if self.identity == None:
            return None
        for i in range(self.order):
            if M[x-ind,i] == self.identity:
                # Right inverse
                return i
            elif M[i,x-ind] == self.identity:
                # Left inverse
                return i
        return None

# Classes for generating specific groups.

class Cyclic_Group(Group):
    """ Returns additive group of integers mod n."""
    def __init__(self, n):
        def op(x,y,p):   return (x+y) % p
        finite_set = [i for i in range(n)]
        # rather than make a new group, we just make this
        # one fit the definition of Group() and super().__init__.
        super().__init__(finite_set, op)

class Multiplicative_Group(Group):
    """ Returns multiplicative group of integers mod n."""
    def __init__(self, n):
        def op(x,y,p):  return (x*y) % p
        self.op = op
        self.finite_set = [i for i in range(1, n)]
        self.order = n-1
        M = self._matrix()
        super().__init__(M)

    def _matrix(self):
        M = np.zeros((self.order, self.order), dtype=int)
        for i, j in enumerate(self.finite_set):
            for m, n in enumerate(self.finite_set):
                M[i,m] = self.op(j, n, self.order+1)
        return M

class Dihedral_Group(Group):
    """ Returns dihedral group of order 2n."""
    def __init__(self, n):
        self.perm_order = n
        self.order = 2*n
        self.op = None
        self.perms = self._gens()
        self.array = self._matrix()
        self.finite_set = allelems(self.array)
        self.index = 0
        self.identity = self._identity()

    def _matrix(self):
        M = np.zeros((self.order, self.order), dtype=int)
        gens = self.perms
        for i, j in enumerate(gens):
            for m, n in enumerate(gens):
                a = np.dot(permtomat(j), permtomat(n))
                a = tuple([list(i).index(1) for i in a])
                M[i,m] = gens.index(a)
        return M

    def _gens(self):
        identity = tuple([i for i in range(self.perm_order)])
        a = [self.perm_order-1] + [i for i in range(self.perm_order-1)]
        b = a[::-1]
        agens, bgens = [tuple(a)], [tuple(b)]
        a, b = permtomat(a), permtomat(b)
        for i in range(self.perm_order-2):
            c = np.dot(permtomat(agens[-1]), a)
            c = tuple([list(i).index(1) for i in c])
            agens.append(c)
        for i in agens:
            c = np.dot(permtomat(i), b)
            c = tuple([list(i).index(1) for i in c])
            bgens.append(c)
        return [identity]+agens+bgens

class Symmetric_Group(Group):
    """ Returns symmetric group of order n!.
        n = 7 has order 5040, expect a long wait.
    """
    def __init__(self, n):
        self.perm_order = n
        self.order = factorial(n)
        self.op = None
        self.perms = list(permutations(range(n)))
        self.array = self._matrix()
        self.finite_set = allelems(self.array)
        self.index = 0
        self.identity = self._identity()
        # No super().__init__ here, we can't reduce this to a
        # set and binary function that Group() can parse.
        # This is due to assigning an arbitrary integer to each
        # permutation for a numerical matrix representation.

    def _matrix(self):
        """ Construct matrix representation of permutations.
            Assigns arbitrary ordering to the set of permutations.
        """
        gens = self.perms
        M = np.zeros((self.order, self.order), dtype=int)
        checked = []
        for i, j in enumerate(gens):
            for m, n in enumerate(gens):
                a = np.dot(permtomat(j), permtomat(n))
                a = tuple([list(i).index(1) for i in a])
                M[i,m] = gens.index(a)
        return M

    def is_abelian(self):
        """ Overrides is_abelian() Group method:
            any symmetric group with n >= 3 is not abelian.
        """
        if self.perm_order >= 3:
            return True
        return False
