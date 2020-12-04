import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from grid_analysis.file_operations import pqr_scraper
from grid_analysis.grid_operations import fit_plane
from grid_analysis.grid_operations import projections

au2volts = 27.211386245988
au2ang = 0.52917721067121

def change_slice_basis(old_bas, xyz_mat):

    print('Cartesian Coordinate Matrix: ', xyz_mat.shape)
    print('Associated Basis:            ', old_bas.shape)

    new_xyz_mat = np.matrix(xyz_mat)*np.matrix(old_bas)

    print('Matrix in the new basis:     ', new_xyz_mat.shape,'\n')

    return new_xyz_mat

def projecter(csv0, basis, mode=None, out=None):

    print('Processing file ---> {0}'.format(csv0),'\n')

    if mode=='esp':
        x, y, z, p = csv0[:,0], csv0[:,1], csv0[:,2], csv0[:,3] * au2volts
        if not out:
            out = "tmp.esp"

        print('array x --->', x)
        print('array y --->', y)
        print('array z --->', z)
        xyz = np.array([x, y, z]).T

    elif mode=='xyz':
        csv0 = np.loadtxt(csv0, skiprows=2)
        x, y, z = csv0[:,1], csv0[:,2], csv0[:,3]
        xyz = np.array([x, y, z]).T

    b1, b2, b3 = basis[0], basis[1], basis[2]
    old_basis = np.array([b1, b2, b3]).T

    new_mat = change_slice_basis(old_basis, xyz)
    new_mat_df = pd.DataFrame(new_mat)

    if mode=='esp' and out:
        new_mat_df[3] = pd.Series(p)
        np.savetxt(out, new_mat_df, fmt='%10.6f %10.6f %10.6f %10.6f')
        return new_mat, out

    elif mode=='xyz' and out:
        np.savetxt(out, new_mat_df, fmt='%10.6f %10.6f %10.6f')
        return new_mat

def change_basis_driver(pqr):
    '''
    Fit plane to a system of pi-atoms (heavy-only) and create
    a new basis of vectors centered on the origin defined as
    the chromophore centroid (all-atoms) projected on the plane.
    The 3 orthonormal basis are the vector normal to the plane,
    a vector belonging to the projections.lane itself and a vector
    normal to this last two.
    '''

    llidx = int(pqr_scraper.find_pqr_linker_lys(pqr))
    x1, y1, z1 = pqr_scraper.extract_ret_pi_pqr(pqr, llidx)
    x2, y2, z2 = pqr_scraper.extract_ret_pqr(pqr, llidx)

    # 1) Plane fit to pi-system
    eq_1 = fit_plane.planefit_3d_points(x1, y1, z1)

    # 2) Whole Retinal Centroid
    ctrd = fit_plane.centroid(x2, y2, z2)

    # 3) Centroid projection onto plane
    ctrd_proj, normal = projections.proj_plane_distance(x1, y1, z1, eq_1[0], eq_1[1], eq_1[2], eq_1[3], ctrd)

    # 4) Atom Choice to select the second rotation axis (C5)
    atom = x1[0], y1[0], y1[0]
    atom_proj = projections.proj_plane_distance(x1, y1, z1, eq_1[0], eq_1[1], eq_1[2], eq_1[3], atom)[0]

    # 5) New Basis
    v1 = normal
    v2 = (atom_proj - ctrd_proj)/np.linalg.norm(atom_proj - ctrd_proj)
    v3 = np.cross(v1, v2)

    return v1, v2, v3, ctrd_proj
