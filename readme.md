# Pyknots

A collection of tools for use with computations in knot theory. Emphasis is on computations involving classical, virtual and singular braids as well as finite racks, quandles, biquandles and more. Pyknots contains information on prime knots up to 13 crossings and connected quandles of order up to 45. These (and any published algorithms involved) are credited in references.txt. Examples of functionality are contained in the examples folder.

## Requirements

Pyknots requires numpy and sympy and is only supported on Python 3.x (Python 2 is untested).

## Setup

Place Pyknots in the desired directory and in Python 3.x do:

>from pyknots import *

This will import the key classes and functions required. Alternatively, relative imports work fine:

>from pyknots.modules.groups import Group

## Introduction

At the core of Pyknots is the Braid class which allows instantiation and manipulation of different types of braids. The files magmas.py, groups.py, and quandles.py contain classes for working with Cayley tables of finite magmas, groups, quandles, biquandles and singquandles, as well as specific objects such as direct products of arbitrary magmas, cyclic groups, dihedral groups, alexander quandles, conjugation quandles and more. notations.py, invariants.py and utils.py contain most of the functions for performing computations on these objects, while structural properties are computed in their respective classes methods. For detailed examples of the following see the examples folder.

### Braids
braids.py contains a class Braid which allows instantiation of a braid object. Input is taken primarily as a list representing the braid notation, however Knots_13.json contains look up tables which allows Braid to identify any string inputs so long as they are in Alexander-Briggs notation, i.e. 3_1 etc. Thanks to the list compiled on M. Saitos webpage (see references.txt).

For more info on notations and knot invariants, see: http://www.indiana.edu/~knotinfo.

### Magmas
magmas.py contains only the class Magma which is a base class for all classes generating Cayley tables. Use this class for instantiating magmas with unknown structure. Magma takes as input a Cayley table in matrix form, which can be a list of lists such as [[a, b], [c, d]] or a numpy array. The input format is not important as all matrices are handled as numpy arrays regardless.

### Groups
groups.py contains the class Group, as well as the Group class constructors: Cyclic_Group, Multiplicative_Group, Symmetric_group, and Dihedral_Group. These inherit from Magma and add methods and attributes specific to group structures. The Group constructor accepts the same input as Magma, or alternatively, a finite set and a well-defined binary group operation (in the form of a python function accepting 2 inputs) from which it defines the corresponding Cayley table itself.

### Quandles
quandles.py contains the class Quandle, Biquandle, and Singquandle, as well as Quandle constructors: Alexander_Quandle, Conj_Quandle, and Trivial_Quandle. Each inherits from Magma class and so accepts a matrix representative of the Cayley table as input. A finite set and quandle operation are not supported, however Pyknots contains data (thanks to L. Vendramin and the GAP package RIG - see references.txt) on the first 630 indecomposable connected quandles with order up to 44 and some of order 45. (The published list contains up to order 47 and Pyknots tables should be updated). To access these, input a string to the quandle constructor indicating the desired quandles number.

### Invariants
invariants.py contains functions for computing specific knot invariants including writhe, alexander polynomials, seifert matrices and more.

### Notations
notations.py contains functions for computing different knot notations for a given braid such as Gauss and Dowker-Thistlethwaite, as well as a more general crossing pattern notation related to the previous two.

### Utils
utils.py contains miscellaneous utility functions for working with permutations, graphs, matrices, and more.
