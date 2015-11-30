#!/usr/bin/env python
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


from glob import glob
import numpy as np, os
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_data import install_data
from Cython.Distutils import build_ext
import subprocess


def runcommand(cmd):
    process = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    c = process.communicate()
    if process.returncode != 0:
        raise Exception("Something went wrong whilst running the command: %s" % cmd)
    return c[0]


class my_install_data(install_data):
    """Add a datadir.txt file that points to the root for the data files. It is
       otherwise impossible to figure out the location of these data files at
       runtime.
    """
    def run(self):
        # Do the normal install_data
        install_data.run(self)
        # Create the file datadir.txt. It's exact content is only known
        # at installation time.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            if '.' not in name:
                destination = os.path.join(libdir, name, "data_dir.txt")
                print "Creating %s" % destination
                if not self.dry_run:
                    f = file(destination, "w")
                    print >> f, self.install_dir
                    f.close()


def find_all_data_files(dn):
    result = []
    for root, dirs, files in os.walk(os.path.join('data', dn)):
        if len(files) > 0:
            files = [os.path.join(root, fn) for fn in files]
            result.append(('share/yaff/' + root[5:], files))
    return result

# Find the MPI arguments required for building the modules.
# Fetch command line, convert to a list, and remove the first item (the command).
# TODO: Make sure this works on most systems...
mpiargs = runcommand('mpicc -show').split()[1:]
# Find arguments needed for FCS
# TODO: find this with pkg-config
fcsargs = "-I/usr/local/include -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lfcs4fortran -lfcs -lfcs_direct -lfcs_ewald -lfcs_fmm -libverbs -lfcs_memd -lfcs_mmm1d -lfcs_mmm2d -lfcs_p2nfft -lgsl -lgslcblas -lfcs_p3m -lfcs_pepc -lfcs_pp3mg -lfcs_vmg -lboost_system -lboost_filesystem -llapack -lcblas -lf77blas -latlas -lfcs_wolf -lfcs_near -lfcs_gridsort -lfcs_resort -lfcs_common -lfcs_pnfft -lfcs_pfft -lfcs_fftw3_mpi -lfcs_fftw3 -lmpichf90 -lgfortran -lgomp -lquadmath -lmpichcxx -lmpich -lopa -lmpl -lrt -lcr -lpthread -lstdc++ -lm -lgcc_s".split()

extralinkargs = mpiargs + fcsargs
extracompileargs = mpiargs + fcsargs

setup(
    name='yaff',
    version='1.0',
    description='YAFF is yet another force-field code.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    package_dir = {'yaff': 'yaff'},
    packages=['yaff', 'yaff/test', 'yaff/pes', 'yaff/pes/test', 'yaff/sampling',
              'yaff/sampling/test', 'yaff/analysis', 'yaff/analysis/test',
              'yaff/tune', 'yaff/tune/test', 'yaff/conversion',
              'yaff/conversion/'],
    cmdclass = {'build_ext': build_ext, 'install_data': my_install_data},
    data_files=[
        ('share/yaff/test', glob('data/test/*.*')),
    ] + find_all_data_files('examples'),
    ext_modules=[
        Extension("yaff.pes.ext",
            sources=['yaff/pes/ext.pyx', 'yaff/pes/nlist.c',
                     'yaff/pes/pair_pot.c', 'yaff/pes/ewald.c',
                     'yaff/pes/dlist.c', 'yaff/pes/grid.c', 'yaff/pes/iclist.c',
                     'yaff/pes/vlist.c', 'yaff/pes/cell.c',
                     'yaff/pes/truncation.c'],
            depends=['yaff/pes/nlist.h', 'yaff/pes/nlist.pxd',
                     'yaff/pes/pair_pot.h', 'yaff/pes/pair_pot.pxd',
                     'yaff/pes/ewald.h', 'yaff/pes/ewald.pxd',
                     'yaff/pes/dlist.h', 'yaff/pes/dlist.pxd',
                     'yaff/pes/grid.h', 'yaff/pes/grid.pxd',
                     'yaff/pes/iclist.h', 'yaff/pes/iclist.pxd',
                     'yaff/pes/vlist.h', 'yaff/pes/vlist.pxd',
                     'yaff/pes/cell.h', 'yaff/pes/cell.pxd',
                     'yaff/pes/truncation.h', 'yaff/pes/truncation.pxd',
                     'yaff/pes/constants.h'],
            include_dirs=[np.get_include()],
            extra_compile_args=extracompileargs,
            extra_link_args=extralinkargs
        ),
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
