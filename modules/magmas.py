""" Contains the base-class Magma, which has all methods
    regarding cayley table properties.

FIXME:
    - fix attributes: decide between elems, finite_set, perms, etc.

"""

import numpy as np
from itertools import permutations, product
from pyknots.modules.utils import issquare, permtomat, allelems, mergedicts, applymorphism
import json

__all__ = ['Magma', 'Direct_Product']

class Magma(object):
    def __init__(self, matrix):
        if issquare(matrix):
            self.array = np.array(matrix)
            self.order = len(matrix[0])
            self.index = np.amin(self.array)
            self.elems = allelems(self.array.tolist())
            self.finite_set = self.elems
        else:
            raise TypeError('Input %s not a square matrix.' % (matrix))

    def __str__(self):
        return str(self.array)

    def pythonic_matrix(self):
        """ Returns 0-indexed matrix """
        vfunc = np.vectorize(lambda x: x-self.index)
        return vfunc(self.array)

    def orbit(self, x, orb=None, checked=None):
        """ The orbit of x in Q is all elements a such that there exists
            b in Q satisfying: x*b == a, as well as the orbit of each a.
        """
        if checked == None:  checked = []
        if orb == None:  orb = []
        checked.append(x)
        ind = self.index
        for i in range(self.order):
            j = self.array[x-ind,i]
            if j not in orb:
                orb.append(j)
        for i in orb:
            if i not in checked:
                orb = self.orbit(i, orb, checked)
        return orb

    def cycle(self, x, cyc=None):
        """ The cycle of an element g is defined as:
            <g> = [g^0, g^1, ... , g^n] where g^n+1 == g
        """
        if cyc == None: cyc = [x]
        ind = self.index
        i = self.array[cyc[-1]-ind, x-ind]
        if i not in cyc:
            cyc.append(i)
            cyc = self.cycle(x, cyc)
        return cyc

    def cayley_digraph(self, gens, color=False):
        """ Returns dict of edges of Cayley digraph. If color=True,
            returns separate dicts for each color arcs, where one color
            corresponds to one generator.
        """
        M, elems = self.array, self.elems
        ind = self.index
        dgraphs = []
        xgens, elem = [gens], []

        for i in gens:
            for j in gens:
                xgens.append([M[i-ind,j-ind]])
        xgens = allelems(xgens)
        for i in xgens:
            elem.append(self.cycle(i))
        elem = allelems(elem)

        for i in gens:
            dgraph = {j: [] for j in elems}
            for j in elem:
                cyc = self.cycle(i, [j])
                if len(cyc) > 1:
                    dgraph[j].append(cyc[1])
            dgraphs.append(dgraph)

        if color:
            return tuple(dgraphs)
        return mergedicts(dgraphs)

    def is_complete(self):
        """ Determines if the matrix will form a complete graph when
            constructing a cayley digraph.
        """
        M = self.array
        elems = self.elems
        reaches = {i: [] for i in elems}
        for i in range(self.order):
            for j in range(self.order):
                reaches[i].append(M[i,j])
        for i in reaches:
            if set(reaches[i]) != set(elems):
                elem = [j for j in elems if j != i]
                if set(reaches[i]) != set(elem):
                    return False
        return True

    def is_left_invertible(self):
        """ X is left invertible if for all a, b in X, there is unique
            x such that x*a = b.
        """
        M = self.array
        ind = self.index
        for i in range(self.order):
            col = []
            for j in range(self.order):
                col.append(M[j,i])
            for c in range(len(col)):
                if not c+ind in col:
                    return False
        return True

    def is_right_invertible(self):
        """ X is right invertible if for all a, b in X, there is unique
            x such that a*x = b.
        """
        # Checks that all elements are in the row.
        M = self.array
        ind = self.index
        for row in M:
            for i in range(self.order):
                if not i+ind in row:
                    return False
        return True

    def is_cyclic(self):
        """ If any element exists such that:
            <g> = [g^0, g^1, ... , g^n] such that g^n+1 == g,
            then the table is cyclic.
        """
        elems = self.elems
        for i in finite_set:
            cyc = self.cycle(i)
            if set(cyc) == set(finite_set):
                return True
        return False

    def is_latin(self):
        """ No element in any row or column is repeated."""
        M = self.array
        for i in range(self.order):
            k = M[i,i]
            for j in range(self.order):
                if i != j:
                    m, n = M[i,j], M[j,i]
                    if n == k or m == k:
                        return False
        return True

    def is_faithful(self):
        """ If any column of the quandle matrix is repeated then the
            table is not faithful.
        """
        M = self.array
        col = []
        for i in range(self.order):
            c = []
            for j in range(self.order):
                c.append(M[j,i])
            if c in col:
                return False
            col.append(c)
        return True

    def is_connected(self):
        """ If one orbit of one element of the table contains all elements
            of the table, then the table is connected.
        """
        M = self.array
        it = np.nditer(M)
        for i in it:
            if len(self.orbit(i)) == self.order:
                return True
        return False

    def is_simple(self):
        raise NotImplemented

    def is_isomorphic(self, other):
        """ Checks all_isomorphisms for positive isomorphic map from
            M1 to M2
        """
        for i in self.all_isomorphisms(other):
            return True
        return False

    def all_isomorphisms(self, other):
        """ Returns iterator of isomorphic maps from M1 to M2 given
            as permutation matrices p such that p**(-1) * M1 * p == M2
        """
        if self.index != other.index:
            M1 = self.pythonic_matrix()
            M2 = other.pythonic_matrix()
            ind = 0
        else:
            M1 = self.array
            M2 = other.array
            ind = self.index
        table = np.arange(self.order)
        sigma = permutations(table)
        for p in sigma:
            M = applymorphism(M1, p)
            if np.array_equal(M, M2):
                yield p

    def all_automorphisms(self):
        for p in self.all_isomorphisms(self):
            yield p

    def automorphism_group(self):
        all_perms = list(self.all_automorphisms())
        """
        all_perms = []
        for i in self.all_isomorphisms(self):
            a = tuple([list(j).index(1) for j in i])
            all_perms.append(a)
        """

        M = np.zeros((len(all_perms), len(all_perms)), dtype=int)
        for i, j in enumerate(all_perms):
            for m, n in enumerate(all_perms):
                a = np.dot(permtomat(j), permtomat(n))
                a = tuple([list(j).index(1) for j in a])
                M[i,m] = all_perms.index(a)
        return M

    def all_inner_automorphisms(self):
        raise NotImplemented

    def is_homomorphic(self, quandle):
        raise NotImplemented

    def all_homomorphisms(self, quandle):
        raise NotImplemented


class Direct_Product(Magma):
    """ Returns table representing the direct product of magma1 and magma2."""
    def __init__(self, magma1, magma2):
        self.magma1, self.magma2 = magma1, magma2
        self.finite_set = list(product(magma1.elems, magma2.elems))
        self.order = magma1.order * magma2.order
        self.elems = [i for i in range(self.order)]
        self.array = self._matrix()
        self.index = 0

    def _matrix(self):
        M1, M2 = self.magma1.array, self.magma2.array
        M = np.zeros((self.order, self.order), dtype=int)
        for i, j in enumerate(self.finite_set):
            for m, n in enumerate(self.finite_set):
                a = (M1[j[0], n[0]], M2[j[1], n[1]])
                M[i,m] = self.finite_set.index(a)
        return M
