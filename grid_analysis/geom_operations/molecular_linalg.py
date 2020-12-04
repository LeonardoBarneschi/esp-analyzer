#!/usr/bin/env python3  
# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function

import os
import scipy
import sys
import time
import numpy as np
import MDAnalysis as mda

from scipy.spatial.distance import cdist
from grid_analysis.file_operations import pqr_scraper
from grid_analysis.geom_operations import vdw_surface
from grid_analysis.geom_operations import getpot
from grid_analysis.geom_operations import getpot_parallel
from grid_analysis.geom_operations import getcoulene_parallel
from grid_analysis.geom_operations import getcoulene

#np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

ang2au = 1/0.52917721067121
au2ang = 0.5291772083

def two_points_dist(a, b):
    c = a - b
    return np.linalg.norm(c)


def dist2sets_cdist(a, b):
    M = cdist(a, b)
    return M

def calc_distance(a,b):
    D=np.empty((a.shape[0],b.shape[0]),dtype=a.dtype)
    for i in nb.prange(a.shape[0]):
        for j in range(b.shape[0]):
            D[i,j]=np.sqrt((a[i,0]-b[j,0])**2+(a[i,1]-b[j,1])**2+(a[i,2]-b[j,2])**2)
    return D

def dist2sets_std(a, b):
    A = np.repeat(a[:,np.newaxis,:], b.shape[0], axis=1)
    B = np.repeat(b[np.newaxis,:,:], a.shape[0], axis=0)
    M = np.sqrt(((A - B)**2).sum(axis=2))
    return M


def pqr2vdwsurf(pqr, density, extrapoints=False, chromophore="RET"):

    idx = pqr_scraper.find_pqr_linker_lys(pqr)
    u = mda.Universe(pqr)

    # TODO GENERALIZE
    sel = u.select_atoms("resname %s or (name CE HE2 HE3 NZ HZ1 and resid %s)" % (chromophore, idx))
    coords = sel.positions
    vdwradii = sel.radii
    elements = sel.types

    surf = vdw_surface.vdw_spheres(coords, elements, 1, density, vdwradii)#, mapper)
    B = np.vstack(surf)
    if extrapoints:
        circlepoints = vdw_surface.multisphereintersec(coords, vdwradii, extrapoints)
        B = np.r_[B, circlepoints]

    D = dist2sets_cdist(coords, B)
    m = D.shape[1]
    F = np.repeat(vdwradii.reshape(-1, 1), m, axis=1)
    tolerance = 0.001
    filtered = D >= F-tolerance
    X = np.all(filtered, axis=0)
    idx = np.where(X == True)
    pts = B[idx[0]]

    return pts

def coulomb_esp(pqr, density=1, scale=1, mode=0, chromophore="RET", surfile=None):

    mode = int(mode)

    idx = pqr_scraper.find_pqr_linker_lys(pqr)
    u = mda.Universe(pqr)
    sel = u.select_atoms("resname %s or (name CE HE2 HE3 NZ HZ1 and resid %s)" % (chromophore, idx))
    vdwradii = sel.radii
    elements = sel.types
    coords = sel.positions

    if not surfile:
        surf = vdw_surface.vdw_spheres(coords, elements, scale, density, vdwradii)
        B = np.vstack(surf)
        D = dist2sets_cdist(coords, B)
        m = D.shape[1]
        F = np.repeat(vdwradii.reshape(-1, 1), m, axis=1)
        tolerance = 0.001
        filtered = D >= F-tolerance
        X = np.all(filtered, axis=0)
        idx = np.where(X == True)
        S = B[idx[0]]

    else:
        S = np.loadtxt(surfile, skiprows=1) * au2ang

    if mode == 0:
        sel = u.select_atoms("all")
        coords = sel.positions
        qp = sel.charges
        #Vp = getpot.getpot(coords, qp, S)
        Vp = getpot_parallel.getpot(coords, qp, S)
        return Vp, S

    # MODIFY !!!
    else:
        qc = sel.charges
        #Vc = getpot.getpot(coords, qc, S)
        Vc = getpot_parallel.getpot(coords, qc, S)
        return Vc, S

def coulinteract(propqr, chrpqr, chromophore="RET"):

    idx = pqr_scraper.find_pqr_linker_lys(propqr)
    pu = mda.Universe(propqr)
    cu = mda.Universe(chrpqr)

    chrsel = cu.select_atoms("resname %s or (name CE HE2 HE3 NZ HZ1 and resid %s)" % (chromophore, idx))
    chrcoords = chrsel.positions
    elements = chrsel.names
    qc = chrsel.charges

    prosel = pu.select_atoms("not (resname %s or (name CE HE2 HE3 NZ HZ1 and resid %s))" % (chromophore, idx))
    procoords = prosel.positions
    qp = prosel.charges
    Vcp = getcoulene_parallel.getcoulene(chrcoords, qc, procoords, qp)

    return Vcp, Vcp.sum(), elements

def write_coulomb_esp(out, V, S, ovw=False):

    if ovw:
        with open(out,'w') as o:
            o.write(str(V.shape[0])+'\n')
            np.savetxt(o, np.c_[S, V], fmt="%10.6f")
    else:
        if not os.path.exists(out):
            with open(out,'w') as o:
                o.write(str(V.shape[0])+'\n')
                np.savetxt(o, np.c_[S, V], fmt="%10.6f")

        else:
            print('\n The file {0} exists. To force overwrite\
                   activate with ovw=True.\n'.format(out))

    return out

if __name__ == "__main__":

    pqr = sys.argv[1]
    out = sys.argv[2]
    mod = sys.argv[3]
    V, XYZ = coulomb_esp(pqr, 5, mod)
    with open(out,'w') as o:
        o.write(str(V.shape[0])+'\n')
        np.savetxt(o, np.c_[XYZ, V], fmt="%10.6f")
