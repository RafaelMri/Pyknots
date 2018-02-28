# Instantiate a braid and perform multiple operations.

from __future__ import print_function
import json
import os
import numpy as np
from sympy import symbols

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
    def __init__(self, braid, scrossings=[], vcrossings=[]):
        b, s, v = self._check(braid, scrossings, vcrossings)
        self.braid, self.scrossings, self.vcrossings = b, s, v
        self.index = max(max(self.braid), abs(min(self.braid)))+1
        self.length = len(self.braid)

    def __str__(self):
        if self.scrossings == [] == self.vcrossings:
            return str(self.braid)
        elif self.scrossings == []:
            return str(self.braid), str(self.vcrossings)
        elif self.vcrossings == []:
            return str(self.braid), str(self.scrossings)
        else:
            return str(self.braid), str(self.scrossings), str(self.vcrossings)

    def __add__(self, B):
        B_braid = [i+self.index for i in B.braid]
        B_scrossings = [i+self.length for i in B.scrossings]
        B_vcrossings = [i+self.length for i in B.vcrossings]
        return Braid(self.braid+B_braid, self.scrossings+B_scrossings, self.vcrossings+B_vcrossings)

    def __mul__(self, B):
        B_scrossings = [i+self.length for i in B.scrossings]
        B_vcrossings = [i+self.length for i in B.vcrossings]
        return Braid(self.braid+B.braid, self.scrossings+B_scrossings, self.vcrossings+B_vcrossings)

    def __eq__(self, B):
        return self.reduce().braid == B.reduce().braid

    #def __ne__(self, B):
    #    return self.reduce().braid != B.reduce().braid

    def __len__(self):
        return self.length

    def _check(self, braid, scrossings, vcrossings):
        if isinstance(braid, list):
            for i in braid:
                if not isinstance(i, int) or i ==0:
                    raise TypeError('Element %s in %s not an integer greater than 0.' %(i, braid))
        elif isinstance(braid, str):
            dirname = os.path.dirname(os.path.realpath(__file__))
            path = os.path.join(os.path.dirname(dirname), 'data', 'knots_13.json')
            with open(path, 'r') as fp:
                knot_dict = json.load(fp)
                try:
                    braid = knot_dict[braid][-1]
                    #self.__init__(braid)
                except KeyError:
                    raise KeyError('Input %s not a valid knot name.' % (braid))
        else:
            raise TypeError('Input %s not type list or string.' % (braid))

        if not isinstance(scrossings, list):
            raise TypeError('Singular crossings %s not a list.' % (scrossings))

        if not isinstance(vcrossings, list):
            raise TypeError('Singular crossings %s not a list.' % (vcrossings))
        return braid, scrossings, vcrossings

    ####################################################
    # braid reduction methods
    ####################################################

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
        b = self.braid
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
        b = self.braid
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

    def _is_free_reduced(self):
        """ For reduce() method, returns True/False depending on whether
            the braid contains any pair [n, -n]. Any braid of this type can
            be reduced from [..., n, n-1, ...] --> [..., ...].
        """
        b = self.braid
        for i, j in enumerate(b):
            if i < len(b)-1 and j + b[i+1] == 0:
                return False
        return True

    def _free_reduce(self):
        """ For reduce() method, removes any pair [n, -n] from the braid word."""
        b = self.braid
        reduced = self._is_free_reduced()
        while not reduced:
            for i, j in enumerate(b):
                if i < len(b)-1 and j + b[i+1] == 0:
                    try:
                        del b[i+1]
                        del b[i]
                    except IndexError:
                        break
            reduced = self._is_free_reduced()
        self.braid = b
        return b

    def is_reduced(self):
        """ For reduce() method, returns True/False depending on whether
            the braid is reduced based on the algorithm given by
            Patrick Dehornoy's paper: "A Fast Method for Comparing Braids."
        """
        b = self.braid
        if not self._is_free_reduced():
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

    def reduce(self):
        """Reduce braid word in place using the combinatorial method outlined
           in Patrick Dehornoy's paper "A Fast Method for Comparing Braids.
           See citations."
        """
        reduced = self.is_reduced()
        while not reduced:
            self.braid = self._free_reduce()
            handle = self._get_handles()
            if handle == []:
                break
            reduced_handle = self._handle_reduce(handle)
            self.braid = self._replace(handle, reduced_handle)
            self.braid = self._free_reduce()
            reduced = self.is_reduced()
        self.__init__(self.braid)
    ######################################################

    def is_singular(self):
        if self.scrossings != []:
            return True
        return False

    def is_virtual(self):
        if self.vcrossings != []:
            return True
        return False

    def is_classical(self):
        if self.scrossings == [] == self.vcrossings:
            return True
        return False

    def start_vars(self, commutes=True):
        """ Generate list of variables to represent starting strands
            of the braid.
        """
        start_vars = []

        for i in range(self.index):
            z = 'z'+str(i)
            z = symbols(z, commutative=commutes)
            start_vars.append(z)

        return start_vars

    def inverse(self):
        """ The braid inverse for any braid b is b' such that
            b*b' = unknot.
        """
        inverse = []
        for i in self.braid[::-1]:
            inverse.append(-i)
        return inverse

    def braid_group(self):
        """ placeholder function. see class docstring."""
        raise NotImplemented

    def apply_func(self, func1, func2=None, vals=None):
        """ Trace braid while applying the given functions at
            each crossing. Function 1 gives x*y. Input function 2 for biquandle
            function.
        """
        if vals is None:
            start_vars = self.start_vars()
        else:
            start_vars = vals
        strands = [i for i in start_vars]

        if func2 is None:
            def func2(x,y): return y

        for i, j in enumerate(self.braid):
            s1, s2 = strands[j-1], strands[j]

            if j > 0:
                strands[j-1], strands[j] = func2(s1, s2), func1(s1, s2)
            elif j < 0:
                strands[j-1], strands[j] = func1(s2, s1), func2(s2, s1)

        return strands

    def apply_singular_func(self, func1, func2, func3, func4=None, vals=None):
        """ Trace singular braid while applying the given functions
            at each crossing. Function 1 gives x*y, function 2 gives
            R_1(x,y), and function 3 gives R_2(x,y). Input function 4 for
            biquandle function.
        """
        if vals is None:
            start_vars = self.start_vars()
        else:
            start_vars = vals
        strands = [i for i in start_vars]

        if func4 is None:
            def func4(x,y): return y

        for i, j in enumerate(self.braid):
            s1, s2 = strands[j-1], strands[j]

            if i in self.scrossings:
                strands[j-1], strands[j] = func2(s1, s2), func3(s1, s2)
            elif j > 0:
                strands[j-1], strands[j] = func4(s1, s2), func1(s1, s2)
            elif j < 0:
                strands[j-1], strands[j] = func1(s2, s1), func4(s2, s1)

        return strands

    def is_colored_by(self, quandle):
        """ Test if the braid is colored by the given Quandle.
            Input is taken as a Quandle, Biquandle, or Singquandle
            class object.

            Temp solution. Should try all possible colorings, since some will
            work while some wont. 1 or more successes means true, no successes 
            means false.
        """
        s = [i for i in range(self.index)]
        if self.is_singular() and quandle.is_singquandle():
            def func1(x,y): return quandle.array1[x,y]
            def func2(x,y): return quandle.array2[x,y]
            def func3(x,y): return quandle.array3[x,y]
            res = self.apply_singular_func(func1, func2, func3, vals=s)

        elif not self.is_singular() and quandle.is_singquandle():
            def func1(x,y): return quandle.array1[x,y]
            res = self.apply_func(func1, vals=s)

        elif not self.is_singular() and quandle.is_biquandle():
            def func1(x,y): return quandle.array1[x,y]
            def func2(x,y): return quandle.array2[x,y]
            res = self.apply_func(func1, func2, vals=s)

        elif not self.is_singular() and quandle.is_quandle():
            def func1(x,y): return quandle.array[x,y]
            res = self.apply_func(func1, vals=s)

        else:
            raise TypeError('Input %s must be a Quandle class object.' %(quandle))

        if res == s:
            return True
        return False
