# YAFF is yet another force-field code
# Copyright (C) 2008 - 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of YAFF.
#
# YAFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# YAFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np

from yaff.topology import Topology
from yaff.ext import Cell


__all__ = ['System']


class System(object):
    def __init__(self, numbers, pos, ffatypes, bonds=None, rvecs=None):
        '''
           **Arguments:**

           numbers
                A numpy array with atomic numbers

           pos
                A numpy array (N,3) with atomic coordinates in bohr.

           ffatypes
                A list of labels of the force field atom types.

           **Optional arguments:**

           bonds
                a numpy array (B,2) with atom indexes (counting starts from
                zero) to define the chemical bonds.

           rvecs
                An array whose rows are the unit cell vectors. At most three
                rows are allowed, each containg three Cartesian coordinates.
        '''
        self.numbers = numbers
        self.pos = pos
        self.ffatypes = ffatypes
        if bonds is None:
            self.topology = None
        else:
            self.topology = Topology(bonds, self.natom)
        self.cell = Cell(rvecs)

    natom = property(lambda self: len(self.pos))
