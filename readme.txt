#Pyknots

A collection of tools for use with computations in knot theory.

###Setup

Place pyknots in your current working directory. Then for quick import, use:

>from pyknots import *

This will import only the key classes and their constructors. Alternatively, relative imports work fine:

>from pyknots.groups import Group

###Introduction

There are 2 (currently) mutually exclusive sides to pyknots, which will be merged at a later date. "braids" provides tools for calculating specific invariants of knots based on their braid notation. "tables", "groups", and "quandles" provide tools for handling Cayley tables using fast numpy arrays. See requirements.txt for info on requirements.

###Braids
braids.py contains a class Braid which allows instantiation of a braid object. Input is taken primarily as
a list representing the braid notation, however Knots_13.json contains look up tables which allows
Braid to identify any string inputs, so long as they are in Alexander-Briggs notation. (ie: 3_1 etc)

For more info on notations and knot invariants, see: http://www.indiana.edu/~knotinfo.
See Braids.py docstring for some examples of the functionality.

###Singbraids
singbraids.py contains the class Singbraid. Input is taken as above, with an optional second input as a list of singular crossing indices (0-indexed!). Then Singbraid([1, 1, 1], [0, 2]) is the braid [1, 1, 1] with the first and last crossing singular. See class docstring for more info.

###Tables
tables.py contains only the class Table, which is a base class for all classes generating Cayley tables. Use this class for instantiating magmas with unknown structure.

###Groups
groups.py contains the class Group, as well as the Group class constructors: Cyclic_Group, Multiplicative_Group, Symmetric_group, and Dihedral_Group. These inherit from Table, and add methods and attributes specific to group structures.

###Quandles
quandles.py contains the class Quandle, Biquandle, and Singquandle, as well as Quandle constructors: Alexander_Quandle, Conj_Quandle, and Trivial_Quandle. Each inherits from Table class as well.

###Utils
utils.py contains miscellaneous utility functions for working with permutations, graphs, matrices, and more. Currently also contains applymorphism, for applying isomorphisms to Cayley tables, as well as directproduct for group products. The location of these functions will most likely change in the future.
