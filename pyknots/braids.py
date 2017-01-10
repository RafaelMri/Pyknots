"""
Instantiate a braid and perform multiple calculations.

EXAMPLES:

Input:
b = Braid([1, -2, 1 -2])
c = Braid('12n_0861')         # loads braid notation from 'knots_13.json' (13 crossings or less)

Methods:
b.alexander_polynomial()      # prints the Alexander polynomial.
b.alexander_matrix()          # prints the associated presentation matrix.
b.dowker_thistlethwaite()     # prints the dowker thistlethwaite notation.
...

Attributes:
b.braid_notation              # prints the braid notation [1, -2, 1, -2]
b.inverse                     # prints the braid inverse
...

"""

from __future__ import print_function
import json
from sympy import Poly
from sympy import primitive
from sympy import symbols
from sympy.abc import x, y, h, t
from sympy.matrices import Matrix
#from sympy.combinatorics.named_groups import SymmetricGroup
from itertools import product
from math import isnan

__all__ = ['Braid']

class Braid(object):
    """Instantiate a braid object by passing Braid() a list of integers
       or a string based on the Knot Table of Invariants. Eg: '12n_712'.
       See citations.

       FIXME:
           - alexander_polynomial() incorrect for '12n_0712' and '12n_0447'

       TODO:
           - More support for links.
           - Possibly separate colorability functions from main class?
           - Alternate notations: Conway etc.
           - Jones polynomial, HOMFLY, Kauffman (can be used to find Jones)
           - Estimate of unkotting number (upper/lower bound)
           - P-colorability returns generator. Is this preferred over a list?
           - permutation_group(): braids are isomorphic to certain
             permutation groups. Potentially inherit sympy SymmetricGroup
             class?
           - braid_group(): B_n is the braid group on n strands. Possibly
             inherit from sympy free group class, or define separate
             braid group class.
    """
    def __init__(self, braid):
        if type(braid) is list:
            if self._check_format(braid):
                self.braid_notation = braid
                self.inverse = self._get_braid_inverse()
                self.braid_index = self._get_braid_index()
                self.braid_length = len(self.braid_notation)
            else:
                raise TypeError('Input %s not type list or string.' % (braid))
        elif type(braid) is str:
            with open('knots_13.json', 'r') as fp:
                knot_dict = json.load(fp)
                try:
                    braid = knot_dict[braid][-1]
                    self.__init__(braid)
                except KeyError:
                    raise TypeError('Input %s not type list or string.' % (braid))

    def __str__(self):
        return str(self.braid_notation)

    def _check_format(self, braid):
        """Check whether braid notation is the proper format,
           will not accept lists containing 0 or non-integers.

           Acceptable examples: [1, -2, 1, -2], [1, 1, 1, 25]
           Not acceptable: [1, 2, a], [0, 12, 3, 10]
        """
        for a in braid:
            if type(a) is int and a != 0:
                continue
            elif a == 0:
                return False
            else:
                return False
        return True

    def _apoly(self, x, y, t, h):
        """ Helper for alexander_polynomial() and others."""
        return t*x + (1-t)*y - h

    def _normalize_apoly(self, apoly):
        """ Helper for alexander_polynomial() method.
            Normalizes alexander polynomial. Factors out greatest
            power of t possible and multiplies by -1 if necessary.
        """
        apoly = primitive(Poly(apoly, t))[-1]
        apoly = apoly.terms_gcd()[-1]
        if apoly.EC() < 0:
            apoly = -1 * apoly
        return apoly.as_expr()

    def _normalize_DT(self, DT):
        """ Helper for dowker_thistlethwaite() method.
            DT notation can be normalized by multiplying -1 to all
            elements. This is because mirror images can not be
            distinguished by DT notation.
        """
        pos = neg = 0
        for i in DT:
            if i > 0:
                pos += 1
            else:
                neg += 1
        if pos >= neg:
            return DT
        else:
            return [-i for i in DT]

    def _get_braid_index(self):
        """ Attribute-setting function.
            Braid index is the total number of strands needed to
            represent the braid.
        """
        self.braid_index = max(max(self.braid_notation),
                               abs(min(self.braid_notation)))+1
        return self.braid_index

    def _get_braid_inverse(self):
        """ Attribute-setting function.
            The braid inverse for any braid b is b' such that
            b*b' = unknot.
        """
        self.inverse = []
        for i in self.braid_notation[::-1]:
            self.inverse.append(-i)
        return self.inverse

    def _get_start_vars(self, commutes=True):
        """ Helper for alexander_system() method. Generate list of
            variables to represent starting strands of the braid.
        """
        start_vars = []
        var_count = 0

        for i in range(self.braid_index):
            # assigns each strand a variable z0 - zn
            z = 'z'+str(i)
            z = symbols(z, commutative=commutes)
            start_vars.append(z)
            var_count += 1

        return start_vars

    def _reduce_system(self, sys, strands, start_vars):
        """ Helper for alexander_system() method."""
        reduced_sys = []

        # closes the braid by mapping the bottom strands to the top
        for equations in sys:
            reduced_sys.append(equations.subs(zip(strands, start_vars)))
        return reduced_sys

    def _is_handle(self, handle):
        """ Helper for reduce() method. Returns a boolean."""
        for i, j in enumerate(handle[1:-1]):
            if abs(j) == abs(handle[0]) or abs(j) == abs(handle[0])-1:
                return False
            for m, n in enumerate(handle[1:-1][::-1]):
                if i == len(handle[1:-1]) - m:
                    break
                elif -j == n:
                    if self._is_handle(handle[i+1:len(handle[1:-1])-m+1]):
                        return False
        return True

    def _replace(self, handle, reduced_handle):
        """ Helper for reduce() method. Replaces handle with reduced
            handle in braid notation b.
        """
        b = self.braid_notation
        for i, j in enumerate(b):
            for m, n in enumerate(b[::-1]):
                if i == len(b) - m:
                    break
                elif -j == n and b[i:len(b)-m] == handle:
                    b = b[:i] + reduced_handle + b[i+len(handle):]
        return b

    def _handle_reduce(self, handle):
        """ Helper for reduce() method. Returns reduced handle. """
        reduced_handle = []
        if handle[0] > 0:
            sign1 = 1
        else:
            sign1 = -1
        for i in handle:
            if abs(i) == abs(handle[0]):
                continue
            elif abs(i) == abs(handle[0])+1:
                if i > 0:
                    sign2 = 1
                else:
                    sign2 = -1
                reduced_handle = reduced_handle + [-sign1*abs(i),
                                                   sign2*(abs(i)-1),
                                                   sign1*abs(i)]
            else:
                reduced_handle.append(i)
        return reduced_handle

    def _get_handles(self):
        """ Helper for reduce() method. Iteratively finds handles in braid. """
        b = self.braid_notation
        handle = []
        for i, j in enumerate(b):
            for m, n in enumerate(b[::-1]):
                if i == len(b) - m:
                    break
                elif -j == n and (len(b)-m)-i >= 3:
                   handle = b[i:len(b)-m]
                   if self._is_handle(handle):
                       return handle
                   else:
                       continue
        return handle

    def det(self, t=-1):
        return self.alexander_polynomial(t=t)

    def crossing_pattern(self):
        """ Returns the crossing pattern defined by tracing the knot from
            a starting point and continuing until the starting point is
            reached again.

            At each crossing we create a list representing
            the under and over arcs. Let 'o' be the overarc, 'u' the
            under arc, then the representation is [u, o, u].

            This representation is used in computing the braids DT code,
            as well as the Gauss code and potentially others.
        """
        strands = [i for i in range(self.braid_index)]
        start_strands = [i for i in strands]
        crossings = []
        valid = True
        while valid:
            for i in self.braid_notation:
                s1, s2 = strands[abs(i)-1], strands[abs(i)]
                if i > 0:
                    strands[i-1], strands[i] = s2, s1
                    crossings.append([s1, s2, s1])
                elif i < 0:
                    strands[abs(i)-1], strands[abs(i)] = s2, s1
                    crossings.append([s2, s1, s2])
            if strands == start_strands:
                valid = False
        return crossings

    def dowker_thistlethwaite(self):
        """ Returns Dowker-Thistlethwaite notation of the braid.

            Dowker Thistlewaite notation, ie DT notation, is a sequence of
            even integers computed by walking the length of the knot from
            a starting point (until the starting point is reached again)
            counting up and labelling each crossing by the current value
            of the count. If the value is even and we are going over a
            crossing, then we take the negative of the value.

            Each crossing will appear exactly twice, and so each is assigned
            a tuple of 1 even and 1 odd integer. The tuples are ordered by
            the odd integers, then the resulting sequence of even is
            the DT notation.
        """
        DT, DT_temp = [], []
        step = 1
        DT_dict = {}
        for i in range(self.braid_length):
            DT_dict[i] = []
        for i, j in enumerate(self.crossing_pattern()):
            crossing = i % self.braid_length
            if 0 == j[0]:
                DT_dict[crossing].append(step)
                step += 1
            elif 0 == j[1]:
                if step % 2 == 0:
                    DT_dict[crossing].append(-step)
                else:
                   DT_dict[crossing].append(step)
                step += 1

        for i in DT_dict:
            temp = []
            for j in DT_dict[i]:
                temp.append(j)
            DT_temp.append(temp)

        for i in range(1, (2*self.braid_length)+1):
            for j, k in enumerate(DT_temp):
                if i % 2 != 0:
                    if i in k:
                        DT.append([m for m in k if m != i])
                    if -i in k:
                        DT.append([m for m in k if m != -i])

        DT = [i for j in DT for i in j]
        DT = self._normalize_DT(DT)
        return DT

    def gauss(self):
        """ Returns Gauss notation of the braid.
        """
        gauss = []
        for i, j in enumerate(self.crossing_pattern()):
            crossing = i % self.braid_length
            crossing += 1
            if 0 == j[0]:
                gauss.append(-crossing)
            elif 0 == j[1]:
                gauss.append(crossing)
        return gauss

    def braid_group(self):
        """ placeholder function. see class docstring."""
        self.braid_group = 'B('+str(self.braid_index)+')'
        #return self.braid_group
        raise NotImplemented

    def permutation_group(self):
        """ placeholder function. see class docstring"""
        #return SymmetricGroup(self.braid_index)
        raise NotImplemented

    def mul(self, B):
        """Multiply braid A by braid B in place.
           Resulting product is equivalent to stacking A over B.
        """
        new = self.braid_notation + B.braid_notation
        self.__init__(new)
        return new

    def invert(self):
        """Call method invert() to invert braid. Inverts braid in place.
           Call attribute inverse to return braid notation of inverse braid."""
        new = self._get_braid_inverse()
        self.__init__(new)
        return new

    def linking_number(self):
        return self.writhe()/2

    def writhe(self):
        """ Find writhe of braid. Writhe = # of positive crossings minus
            # of negative crossings.
        """
        writhe = 0
        for i in self.braid_notation:
            if i < 0:
                writhe -= 1
            else:
                writhe += 1
        return writhe

    def is_free_reduced(self):
        """ For reduce() method, returns True/False depending on whether
            the braid contains any pair [n, -n]. Any braid of this type can
            be reduced from [..., n, n-1, ...] --> [..., ...].
        """
        b = self.braid_notation
        for i, j in enumerate(b):
            if i < len(b)-1 and j + b[i+1] == 0:
                return False
        return True

    def is_reduced(self):
        """ For reduce() method, returns True/False depending on whether
            the braid is reduced based on the algorithm given by
            Patrick Dehornoy's paper: "A Fast Method for Comparing Braids."
        """
        b = self.braid_notation
        if not self.is_free_reduced():
            return False
        for i, j in enumerate(b):
            for m, n in enumerate(b[::-1]):
                if i == len(b) - m:
                    break
                elif -j == n and (len(b)-m)-i >= 3:
                   handle = b[i:len(b)-m]
                   if self._is_handle(handle):
                       return False
        return True

    def free_reduce(self):
        """ For reduce() method, removes any pair [n, -n] from the braid word."""
        b = self.braid_notation
        reduced = self.is_free_reduced()
        while not reduced:
            for i, j in enumerate(b):
                if i < len(b)-1 and j + b[i+1] == 0:
                    try:
                        del b[i+1]
                        del b[i]
                    except IndexError:
                        break
            reduced = self.is_free_reduced()
        self.braid_notation = b
        return b

    def reduce(self):
        """Reduce braid word in place using the combinatorial method outlined
           in Patrick Dehornoy's paper "A Fast Method for Comparing Braids.
           See citations."
        """
        reduced = self.is_reduced()
        while not reduced:
            self.braid_notation = self.free_reduce()
            handle = self._get_handles()
            if handle == []:
                break
            reduced_handle = self._handle_reduce(handle)
            self.braid_notation = self._replace(handle, reduced_handle)
            self.braid_notation = self.free_reduce()
            reduced = self.is_reduced()
        self.__init__(self.braid_notation)
        return self.braid_notation

    def _apoly_kwargs_wrapper(self, kwargs):
        """ Alexander polynomial machinery kwarg wrapper. """
        if kwargs == {}:
            return None, None
        else:
            if 't' in kwargs:
                t_value = kwargs['t']
            else:
                t_value = None
            if 'modulus' in kwargs:
                mod = kwargs['modulus']
            else:
                mod = None
        return t_value, mod

    def alexander_system(self, **kwargs):
        """ Create a system of equations for each crossing of the braid.
            Accepts keyword arguments t.
        """
        t_value, mod = self._apoly_kwargs_wrapper(kwargs)

        sys, all_vars = [], []
        start_vars = self._get_start_vars()
        strands = [i for i in start_vars]
        var_count = len(start_vars)

        for i in self.braid_notation:
            z = 'z'+str(var_count)
            z = symbols(z)
            all_vars.append(z)
            var_count += 1
            if i > 0:
                s1, s2 = strands[i-1], strands[i]
                strands[i], strands[i-1] = z, s2
                if not t_value:
                    sys.append(self._apoly(s1, s2, t, z))
                else:
                    sys.append(self._apoly(s1, s2, t_value, z))
            elif i < 0:
                i = abs(i)
                s1, s2 = strands[i-1], strands[i]
                strands[i-1], strands[i] = z, s1
                if not t_value:
                    sys.append(self._apoly(z, s1, t, s2))
                else:
                    sys.append(self._apoly(z, s1, t_value, s2))

        all_vars = start_vars + all_vars
        self._all_vars = [a for a in all_vars if not a in strands]
        sys = self._reduce_system(sys, strands, start_vars)
        return sys

    def alexander_matrix(self, **kwargs):
        """ Returns presentation matrix of alexander polynomial.
            Accepts keyword arguments t, modulus.
        """
        # Specifying t and modulus not working. See class docstring for info
        t_value, mod = self._apoly_kwargs_wrapper(kwargs)
        sys = self.alexander_system()
        coeffs = []

        for i in sys:
            row = []
            for var in self._all_vars:
                if not t_value:
                    row.append(i.coeff(var))
                if t_value is not None:
                    c = Poly(i.coeff(var), t).subs(t, t_value)
                    if not mod:
                        row.append(c.as_expr())
                    else:
                        c = c.as_expr() % mod
                        row.append(c)
            coeffs.append(row)

        alexander_matrix = Matrix(coeffs)
        return alexander_matrix

    def alexander_polynomial(self, **kwargs):
        """ Calculate Alexander Polynomial. Accepts keyword
            arguments t, modulus.
        """
        t_value, mod = self._apoly_kwargs_wrapper(kwargs)

        alexander_matrix = self.alexander_matrix()
        matrix_copy = alexander_matrix[:,:]
        matrix_copy.col_del(-1)
        matrix_copy.row_del(-1)
        apoly = matrix_copy.det()

        # Sympy det() method sometimes returns 0 or NaN on large
        # polynomials,the following are some attempts to weed these out.
        # If we run reduce() on the braid and then retry to calculate the
        # alexander polynomial, it often fixes the problem. See the
        # class docstring for more info on known problems.

        try:
            if isnan(apoly) or apoly == 0:
                if not self.is_reduced():
                    self.reduce()
                    self.alexander_polynomial()
        except:
            pass

        apoly = self._normalize_apoly(apoly)
        if not t_value and mod is not None:
            coeffs = Poly(apoly).all_coeffs()
            coeffs = [i % mod for i in coeffs]
            apoly = Poly.from_list(coeffs, t).as_expr()
        elif t_value is not None and mod is not None:
            apoly = apoly.subs(t, t_value) % mod
        elif t_value is not None and not mod:
            apoly = apoly.subs(t, t_value)
        return apoly

    def p_colorable(self, p=50, *args, **kwargs):
        """ Returns generator of alexander quandles that will color
            the braid. Working on better implementation. See class
            docstring.
        """
        apoly = self.alexander_polynomial()
        for i in range(2, p):
            for j in range(1, i):
                if apoly.subs(t, j) % i == 0:
                    yield (i, j)

    # Working on better implementation. See class docstring.
    #
    #def p_colorable(self, p=50, *args, **kwargs):
    #    """ Returns alexander quandles that will color the braid.
    #        If no value for p is specified, 50 is used.
    #        Output: (p, a) for Z_p[t]/(t - a).
    #    """
    #
    #    gen = self._p_colorable(p, *args, **kwargs)
    #    return [i for i in gen]

    def _coloring_kwargs_wrapper_(self, kwargs):
        """ Braid coloring machinery kwarg wrapper. """
        if kwargs == {}:
            t = -1
            p = self.det(t)
        else:
            if 't' in kwargs:
                t = kwargs['t']
            else:
                t = -1
            if 'p' in kwargs:
                p = kwargs['p']
            else:
                p = self.det()
        return t, p

    def coloring_system(self, **kwargs):
        """ Returns the system generated by coloring the strands
            by the Fox coloring. Accepts keyword arguments t, p,
            to color by different quandles.
        """
        t, p = self._coloring_kwargs_wrapper_(kwargs)
        all_vars = self._get_start_vars()
        strands = [i for i in all_vars]

        for i in self.braid_notation:
            if i > 0:
                s1, s2 = strands[i-1], strands[i]
                strands[i-1], strands[i] = s2, self._apoly(s1, s2, t, 0)
            elif i < 0:
                i = abs(i)
                s1, s2 = strands[i-1], strands[i]
                strands[i-1], strands[i] = self._apoly(s2, s1, t, 0), s1

        for i, j in enumerate(strands):
            strands[i] = j - all_vars[i]

        return strands

    def coloring_matrix(self, **kwargs):
        """ Returns the matrix generated by coloring the strands
            by the Fox coloring. Accepts keyword arguments t, p,
            to color by different quandles.
        """
        t, p = self._coloring_kwargs_wrapper_(kwargs)
        sys = self.coloring_system(t=t, p=p)
        start_vars = self._get_start_vars()
        coeffs = []

        for equations in sys:
            row = []
            for var in start_vars:
                coef = equations.coeff(var) % p
                row.append(coef)
            coeffs.append(row)

        M = Matrix(coeffs)
        return M

    def coloring_basis(self, **kwargs):
        """ Returns the basis generated by coloring the strands
            by the Fox coloring. The set of colors formed by linear
            combinations of this basis will generate all possible
            colorings of the braid. Accepts keyword arguments t, p,
            to color by different quandles.
        """
        t, p = self._coloring_kwargs_wrapper_(kwargs)

        M = self.coloring_matrix(t=t, p=p)
        basis = []

       # Maybe one day sympy will implement matrix operations on
       # finite fields? Until then we will settle with trying all
       # possibilities.
       #
       # Note: Sage implements the Block-Wiedemann algorithm for computing
       # the nullspace on a finite field using Cython C extensions and the
       # C library Givaro.

        for i in product(range(p), repeat=self.braid_index-2):
            k = i + (0, 1)
            b = M.multiply(Matrix(k)).applyfunc(lambda x: x % p)
            if b.is_zero:
                basis.append(Matrix(k))
                basis.append(Matrix(k).applyfunc(lambda x: (1-x)%p))
                return basis
        return basis
