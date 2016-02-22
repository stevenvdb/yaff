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

import numpy as np

from yaff.log import log, timer
from yaff.pes.splines import PairPotSpline
from yaff.pes.ff import ForcePart, ForcePartPair, ForceField
from molmod.units import angstrom, kjmol

__all__ = [
    'swap_pairpot_splines','ForcePartSpline'
]

def swap_pairpot_splines(ff, **kwargs):
    '''
    Replace all the ForcePartPair parts of a force field with one tabulated
    ForcePartSpline potential.
    '''
    #TODO Does not work well when no truncation is applied...
    #TODO Checking the scalings is ugly...
    newff = ForceField(ff.system, [], nlist=ff.nlist)
    scaling_groups = []
    part_lists = []
    for part in ff.parts:
        if isinstance(part, ForcePartPair):
            scales = np.array([part.scalings.scale1,part.scalings.scale2,part.scalings.scale3])
            # Check if we encountered a part with exactly the same scaling rules before...
            found = False
            for igroup, scale_group in enumerate(scaling_groups):
                #...yes! Add this part to the existing list of parts with these scaling rules
                if np.all(scales==scale_group):
                    part_lists[igroup].append(part)
                    found = True
            #...no! Create a new list that will contain all parts with these scaling rules
            if not found:
                scaling_groups.append(scales)
                part_lists.append([part])
        else:
            newff.add_part(part)
    for i in xrange(len(part_lists)):
        pairpotff = ForceField(ff.system, [], nlist=ff.nlist)
        for part in part_lists[i]:
            pairpotff.add_part(part)
        pairpotff.compute()
        pairpotspline = ForcePartSpline(pairpotff, suffix='_%d'%i, **kwargs)
        newff.add_part(pairpotspline)
    return newff

def get_scaling(neigh,stab):
    if neigh['r0'] == 0 and neigh['r1'] == 0 and neigh['r2'] == 0:
        a = neigh['a']
        b = neigh['b']
        for scale in stab:
            if scale['a'] == a:
                if scale['b'] == b:
                    return scale['scale']
            if scale['a'] == b:
                if scale['b'] == a:
                    return scale['scale']
    return 1.0

class ForcePartSpline(ForcePart):
    '''A pairwise (short-range) non-bonding interaction term calculated by
       spline interpolation from tabulated values.
    '''
    def __init__(self, ncff, rmin=1.0*angstrom, rmax=None, npoints=1000, suffix='', store=None, load=None):
        ForcePart.__init__(self, 'splines'+suffix, ncff.system)
        assert ncff.nlist is not None
        assert not ((load is not None) and (store is not None))
        self.system = ncff.system
        self.nlist = ncff.nlist
        self.scalings = ncff.parts[0].scalings
        scales = np.array([self.scalings.scale1,self.scalings.scale2,self.scalings.scale3])
        for part in ncff.parts:
            assert isinstance(part, ForcePartPair), "The spline potential can \
                only be used for short-ranged non-bonded interactions"
            this_scales = np.array([part.scalings.scale1,part.scalings.scale2,part.scalings.scale3])
            assert np.all(scales==this_scales), "All ForcePartPair need to have \
                the same scalings in order to be represented by one spline potential"
        if rmax is None: rmax = self.nlist.rcut*1.05
        assert rmax >= self.nlist.rcut
        with log.section('SPLINES'):
            log('Constructing splines for %i parts:&%s.' % (
                len(ncff.parts), ', '.join(part.name for part in ncff.parts)
            ))
        if load is not None:
            phi = np.loadtxt('%s_phi%s.dat'%(load,suffix))
            dphi = np.loadtxt('%s_dphi%s.dat'%(load,suffix))
            nsplines = phi.shape[0]
            x = np.linspace(rmin, rmax, num=1000, endpoint=True)
            nffa = np.amax(ncff.system.ffatype_ids)+1
        else:
            import matplotlib.pyplot as plt
            with timer.section('MAKE %s' % self.name):
                # Construct x-y pairs
                x = np.linspace(rmin, rmax, num=1000, endpoint=True)
                y = np.zeros(x.shape)
                dy = np.zeros(x.shape)
                ncff.nlist.nneigh = 1
                nffa = np.amax(ncff.system.ffatype_ids)+1
                nsplines = (nffa*(nffa+1))/2
                stab = ncff.parts[0].scalings.stab.copy()
                phi = []
                dphi = []
                gpos = np.zeros(self.system.pos.shape)
                for iffa in xrange(nffa):
                    for jffa in xrange(iffa,nffa):
                        a = np.where(ncff.system.ffatype_ids==iffa)[0][0]
                        # TODO Check that scale of interaction between a and b is 1
                        b = np.where(ncff.system.ffatype_ids==jffa)[0][-1]
                        #if a<b: b,a = a,b
                        for i in xrange(x.shape[0]):
                            gpos[:] = 0.0
                            ncff.nlist.neighs[0] = (a, b, x[i], 0.0, 0.0, x[i], 1, 0, 0)
                            y[i] = ncff.compute(gpos=gpos)
                            dy[i] = -gpos[a,2]
                        phi.append(y.copy())
                        dphi.append(dy.copy())
                        index = nffa*iffa - (iffa*(iffa-1))/2 + jffa - iffa
                        name = "%s-%s" % (ncff.system.ffatypes[iffa],ncff.system.ffatypes[jffa])
                        scale = get_scaling(ncff.nlist.neighs[0],stab)
                        #assert scale==1.0
                        #print index, name, a, b, np.amax(np.abs(y))/kjmol
                        #plt.clf()
                        #plt.plot(x/angstrom, y/kjmol)
                        #plt.savefig("%05d-%s.png" % (index,name))
                phi = np.asarray(phi)
                dphi = np.asarray(dphi)
            if store is not None:
                np.savetxt('%s_phi%s.dat'%(store,suffix),phi)
                np.savetxt('%s_dphi%s.dat'%(store,suffix),dphi)
                np.savetxt('%s_x%s.dat'%(store,suffix),x)
        self.spline_pot = PairPotSpline(nsplines, self.system.ffatype_ids, nffa, x, phi, dphi, self.scalings.stab)
        self.nlist.update()

    def _internal_compute(self, gpos, vtens, hess):
        assert hess is None
        with timer.section('PP %s' % self.name):
            return self.spline_pot.compute(self.nlist.neighs, gpos, vtens, self.nlist.nneigh)
