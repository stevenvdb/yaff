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


import shutil, tempfile

from yaff import *



def check_consistent(pf1, pf2):
    assert len(pf1.sections) == len(pf2.sections)
    for prefix1, section1 in pf1.sections.iteritems():
        section2 = pf2[prefix1]
        assert section1.prefix == section2.prefix
        assert len(section1.definitions) == len(section2.definitions)
        for suffix1, definition1 in section1.definitions.iteritems():
            definition2 = section2.definitions[suffix1]
            assert len(definition1.lines) == len(definition2.lines)
            for (counter1, data1), (counter2, data2) in zip(definition1.lines, definition2.lines):
                assert data1 == data2


def test_consistency_io():
    for fn_parameters in 'test/parameters_bks.txt', 'test/parameters_water.txt':
        pf1 = Parameters.from_file(context.get_fn(fn_parameters))
        dirname = tempfile.mkdtemp('yaff', 'test_consistency_parameters')
        try:
            pf1.write_to_file('%s/parameters_foo.txt' % dirname)
            pf2 = Parameters.from_file('%s/parameters_foo.txt' % dirname)
            check_consistent(pf1, pf2)
        finally:
            shutil.rmtree(dirname)


def test_consistency_copy():
    for fn_parameters in 'test/parameters_bks.txt', 'test/parameters_water.txt':
        pf1 = Parameters.from_file(context.get_fn(fn_parameters))
        pf2 = pf1.copy()
        check_consistent(pf1, pf2)


def test_from_file_bks():
    fn_pars = context.get_fn('test/parameters_bks.txt')
    pf = Parameters.from_file(fn_pars)
    assert pf['EXPREP'].complain.filename == fn_pars
    assert pf['EXPREP']['CPARS'].complain.filename == fn_pars
    assert pf['EXPREP']['CPARS'][0][1] == '       O        O  1.3887730000e+03  2.7600000000e+00'
    assert pf['DAMPDISP']['UNIT'][2][0] == 10
    assert pf['FIXQ']['SCALE'][-1][1] == '3 1.0'
    assert len(pf['FOO'].definitions) == 0
    assert len(pf['FOO']['BARR'].lines) == 0


def test_from_file_water_2():
    fn_pars1 = context.get_fn('test/parameters_water_dampdisp1.txt')
    fn_pars2 = context.get_fn('test/parameters_water_exprep1.txt')
    pf = Parameters.from_file([fn_pars1, fn_pars2])
    assert pf['DAMPDISP'].complain.filename == fn_pars1
    assert pf['DAMPDISP']['SCALE'].complain.filename == fn_pars1
    assert pf['DAMPDISP']['SCALE'][0][0] == 4
    assert pf['DAMPDISP']['SCALE'][0][1] == '1 1.0'
    assert pf['EXPREP'].complain.filename == fn_pars2
    assert pf['EXPREP']['PARS'].complain.filename == fn_pars2
    assert pf['EXPREP']['PARS'][0][0] == 8
    assert pf['EXPREP']['PARS'][0][1] == 'O 4.2117588157e+02 4.4661933834e+00'


def test_complain():
    complain = Complain('foo.bar')
    try:
        complain(22, 'Warning! Warning!')
        assert False
    except IOError:
        pass
    try:
        complain(None, 'High voltage!')
        assert False
    except IOError:
        pass
