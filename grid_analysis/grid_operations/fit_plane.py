#!/usr/bin/env python3 

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from grid_analysis.file_operations import pqr_scraper

au2ang = 0.5291772083
ang2au = 1/au2ang

def planefit_3d_points(x, y, z):
    '''
    Ordinary least squares regression
    '''

    tmp_A = []
    tmp_b = []

    for i in range(len(x)):
        tmp_A.append([x[i], y[i], 1])
        tmp_b.append(z[i])

    b = np.matrix(tmp_b).T
    a = np.matrix(tmp_A)
    fit = (a.T * a).I * a.T * b
    errors = b - a * fit
    residual = np.linalg.norm(errors)

    A, B, C, D = (fit[0]/-fit[2])[0, 0],(fit[1]/-fit[2])[0, 0],(-1/-fit[2])[0, 0],-1

    return A, B, C, D

def norm2(X):
	return np.sqrt(np.sum(X ** 2))

def normalized(X):
	return X / norm2(X)

def point_dist(p1, p2):

    p3 = p1 - p2

    #print(np.linalg.norm(p3))

    return np.linalg.norm(p3)

def centroid(x, y, z):
    xmax, xmin = max(x), min(x)
    xmid = (xmax + xmin)/2
    ymax, ymin = max(y), min(y)
    ymid = (ymax + ymin)/2
    zmax, zmin  = max(z), min(z)
    zmid = (zmax + zmin)/2
    centroid = np.asarray([xmid, ymid, zmid])

    return centroid

def plot_planefit(x, y, z, A, B, C, D):

    print(A, B, C, D)
    # plot raw data
    plt.figure()
    ax = plt.subplot(111, projection='3d')
    ax.scatter(x, y, z, color='b')

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    X,Y = np.meshgrid(np.linspace(xlim[0], xlim[1], 40),
                      np.linspace(ylim[0], ylim[1], 40))

    Z = np.zeros(X.shape)

    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r,c] = (A * X[r,c] + B * Y[r,c] + D)/-C

    ax.scatter(X,Y,Z)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show()

def create_grid(origin, a1, a2, a3, xd, yd, zd, out):

    # 84 / 113 / 132 / 161 / 192 
    nbin = 84
    print('Check for basis orthogonality')
    print('{:.1f}'.format(np.dot(a1, a2)))
    print('{:.1f}'.format(np.dot(a1, a3)))
    print('{:.1f}'.format(np.dot(a2, a3)), '\n')

    grid = []

    # no need to change basis
    for i in np.linspace(-xd/2, xd/2, nbin):
        for j in np.linspace(-yd/2, yd/2, nbin):
            for k in np.linspace(-zd/2, zd/2, nbin):
                grid.append(origin + i*a1 + j*a2 + k*a3)

    grid = np.array(grid)

    with open(out, 'w') as f:
        f.write("%d \n\n" % (len(grid)))
        np.savetxt(f, np.c_[np.ones(len(grid)), grid], fmt='%3d %10.6f %10.6f %10.6f')

    return

def slice_grid(origin, a1, a2, a3, xd, yd, zd, out, wfn):

    nbin = 84
    slice_= []

    # no need to change basis
    for j in np.linspace(-yd/2, yd/2, nbin):
        for k in np.linspace(-zd/2, zd/2, nbin):
            slice_.append(origin + j*a2 + k*a3)

    slice_ = np.array(slice_)

    if out:
        with open(out,'w') as f:
            print(out)
            f.write("%d \n\n" % (len(slice_)))
            np.savetxt(f, np.c_[np.ones(len(slice_)), slice_], fmt='%3d %10.6f %10.6f %10.6f')

    if wfn:
        with open(wfn,'w') as f:
            f.write("%d \n\n" % (len(slice_)))
            np.savetxt(f, np.c_[slice_ * ang2au], fmt='%10.6f %10.6f %10.6f')

    return np.c_[slice_]

def fit_piplane_driver(pqr):

    llidx = int(pqr_scraper.find_pqr_linker_lys(pqr))
    print(llidx)
    x, y, z = pqr_scraper.extract_ret_pi_pqr(pqr, llidx)
    #plane_eq = fit_plane.planefit_3d_points(x, y, z)
    plane_eq = planefit_3d_points(x, y, z)
    return plane_eq

if __name__ == "__main__":

    pass
