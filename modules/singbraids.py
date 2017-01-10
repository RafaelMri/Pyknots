from pyknots.modules.braids import Braid
from sympy import symbols

__all__ = ['Singbraid']

class Singbraid(object):
    """ Instantiate singular braid object by passing Singbraid() a
        list of integers, as well as a second list of which crossings
        are singular (0-indexed). Default behavior is to treat all
        crossings as classical.
    """
    def __init__(self, braid, singcrossings=None):
        if singcrossings is None:
            #if isinstance(braid, list):
            #    singcrossings = [i for i in range(len(braid))]
            singcrossings = []
        if self._check(braid, singcrossings):
            self.braid = Braid(braid)
            self.braid_notation = braid
            self.singcrossings = singcrossings
        else:
            raise TypeError('Input %s, %s not type list.' % (matrix, singcrossings))

    def _check(self, braid, singcrossings):
        if not isinstance(braid, list):
            return False
        if not isinstance(singcrossings, list):
            return False
        return True

    #def _solve(self, strands, start_vars):
    #    sys = []
    #    for i, j in enumerate(strands):
    #        expr = start_vars[i] - j
    #        sys.append(expr)
    #    return sys


    #def is_singular(self):
    #    """ Check if the braid has a singular crossing."""
    #    if self.singcrossings:
    #        return True
    #    return False

    def is_singular(self, x):
        """ Check whether x is a singular crossing in the braid."""
        if x in self.singcrossings:
            return True
        return False

    def weight(self):
        """ Computes weight (related to linking number). -- experimental."""
        weight = 0
        for i, j in enumerate(self.braid_notation):
            if not self.is_singular(i):
                if j < 0:
                    weight += 1
                else:
                    weight -= 1
        return weight

    def alex_singquandle(self, a=-1, modulus=None):
        """ Return alexander singular quandle."""
        strands = self.alex_coloring(a=a)
        start_vars = self.braid._get_start_vars()
        M = []
        for i in strands:
            row = []
            for j in start_vars:
                row.append(i.coeff(j))
            M.append(row)
        return M

    def free_group(self, word1, word2, word3):
        """ Input words of Sympy symbols and compute free group relations."""
        def func1(x,y): return word1
        def func2(x,y): return word2
        def func3(x,y): return word3
        return self.singular_applyfunc(func1, func2, func3, commutes=False)

    def core_conjugation(self):
        """ Returns free group relations given by core conjugation."""
        def func1(x,y): return y*(x**(-1))*y
        def func2(x,y): return x*(y**(-1))*x
        def func3(x,y): return x
        return self.singular_applyfunc(func1, func2, func3, commutes=False)

    def alex_coloring(self, a=-1):
        """ Returns relations given by alexander coloring of singular braids."""
        def func1(x,y): return (2*y - x)
        def func2(x,y): return (a*x + (1-a)*y)
        def func3(x,y): return ((a-1)*x + (2-a)*y)
        return self.singular_applyfunc(func1, func2, func3, commutes=True)

    def singular_applyfunc(self, func1, func2, func3, commutes=True):
        """ Apply arbitrary functions to braid crossings and return end result.
            Use commutes=False for free group operations.
        """
        B = self.braid_notation
        start_vars = self.braid._get_start_vars(commutes=commutes)
        all_vars = start_vars.copy()
        strands = start_vars.copy()
        var_count = len(start_vars)

        for i, j in enumerate(B):
            z = 'z'+str(var_count)
            z = symbols(z, commutative=False)
            all_vars.append(z)
            var_count += 1

            s1, s2 = strands[j-1], strands[j]
            if self.is_singular(i):
                strands[j-1], strands[j] = func2(s1, s2), func3(s1, s2)
            elif j > 0:
                strands[j-1], strands[j] = s2, func1(s1, s2)
            elif j < 0:
                strands[j-1], strands[j] = func1(s2, s1), s1

        self._all_vars = all_vars
        #res = self._solve(strands, start_vars)
        return strands


# TESTS:
#s = Singbraid([1, 1]) # Hopf Link. Example 4.4 in [1]
