# -*- coding: utf-8 -*-
# YAFF is yet another force-field code
# Copyright (C) 2011 - 2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
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
#--
'''Defines the context in which Yaff is used.

   This module controls global parameters that are purely technical, e.g. the
   location of data files. It is certainly not meant to keep track of input
   parameters for a computation.

   This module contains a context object, an instance of the :class:`Context`
   class. For now, its functionality is rather limited. It tries to figure
   out the location of the data directory. If it is not specified in the
   environment variable ``YAFFDATA``, it is assumed that the data is located
   in a directory called ``data``. If the data directory does not exist, an
   error is raised.
'''


import os
from glob import glob


__all__ = ['context', 'Context']


class Context(object):
    '''Finds out where the data directory is located etc.

       The data directory contains data files with standard basis sets and
       pseudo potentials.
    '''
    def __init__(self):
        self.data_dir = os.getenv('YAFFDATA')
        if self.data_dir is None:
            fn_data_dir = os.path.join(os.path.dirname(__file__), 'data_dir.txt')
            if os.path.isfile(fn_data_dir):
                with open(fn_data_dir) as f:
                    self.data_dir = os.path.join(f.read().strip(), 'share/yaff')
        if self.data_dir is None:
            self.data_dir = './data'
        self.data_dir = os.path.abspath(self.data_dir)
        if not os.path.isdir(self.data_dir):
            raise IOError('Can not find the data files. The directory %s does not exist.' % self.data_dir)

    def get_fn(self, filename):
        '''Return the full path to the given filename in the data directory.'''
        return os.path.join(self.data_dir, filename)

    def glob(self, pattern):
        '''Return all files in the data directory that match the given pattern.'''
        return glob(self.get_fn(pattern))


context = Context()
