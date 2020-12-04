from __future__ import division, absolute_import, print_function

import math
import numpy as np
import sys
import time
import itertools

from grid_analysis.file_operations import pqr_scraper
from grid_analysis.geom_operations import molecular_linalg
from grid_analysis.file_operations import file_io

"""
A module to generate van der Waals surfaces of molecules.
"""

def spheres_intersec(o1, o2, r1, r2, pts):

    # plane equation (intersecion)
    A, B, C, D = 2*(o2[0]-o1[0]),\
                 2*(o2[1]-o1[1]),\
                 2*(o2[2]-o1[2]),\
                 o1[0]**2 - o2[0]**2 + o1[1]**2 - o2[1]**2 + o1[2]**2 - o2[2]**2 - r1**2 + r2**2

    # normal to the plane                              
    n = np.array([A, B, C])

    d = molecular_linalg.two_points_dist(o1, o2)

    # equation for intersection radius (it is a circle)
    alfa = np.arccos( ( (r1**2 + d**2 - r2**2) / (2*r1*d) ) )

    #cutoff = 0.075
    cutoff = 0.05

    # line - plane intersection
    t = (A*o1[0] + B*o1[1] + C*o1[2] + D) / (A*(o1[0]-o2[0]) + B*(o1[1]-o2[1]) + C*(o1[2]-o2[2]))
    x = o1[0] + t*(o2[0] - o1[0])
    y = o1[1] + t*(o2[1] - o1[1])
    z = o1[2] + t*(o2[2] - o1[2])

    o3 = x, y, z

    r = r1*np.sin(alfa) + cutoff

    px, py = 1, 1
    pz = (A*px + B*py + D)/(-C)
    cx, cy = 2, -1
    cz = (A*cx + B*cy + D)/(-C)

    ref1 = np.array([px, py, pz]) - np.array([o3[0], o3[1], o3[2]])
    v1 = ref1 / np.linalg.norm(ref1)
    v2 = np.cross(n, v1)
    v2 = v2 / np.linalg.norm(v2)

    theta = np.linspace(0, 2*np.pi, pts)
    circle = []
    for phi in theta:
        newc = o3 + r*np.cos(phi)*v1 + r*np.sin(phi)*v2
        circle.append(newc)

    return np.array([circle]).reshape(-1, 3)

def multisphereintersec(nucleii, vdwradii, pts):

    A = np.c_[nucleii, vdwradii]
    for j in itertools.combinations(A, 2):
        arr1, arr2 = j
        o1r1 = [arr1[:-1], arr1[-1]]
        o2r2 = [arr2[:-1], arr2[-1]]
        circle = spheres_intersec(o1r1[0], o2r2[0], o1r1[1], o2r2[1], pts)
        try:
            circles = np.r_[circles, circle]
        except:
            circles = np.array(circle)

    return circles

def fibonacci_sphere(s=1):

    # Leo Vectorization
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
    p = np.full((s, ), phi)
    j = np.full((s, ), s-1)
    s = np.arange(0, s, 1)
    y = 1-(s/j)*2
    x = np.cos(p*s) * np.sqrt(1 - y*y)
    z = np.sin(p*s) * np.sqrt(1 - y*y)
    M = np.array([x, y, z]).T
    return M

def surface(n):
    """Computes approximately n points on unit sphere. Code adapted from GAMESS.
    Parameters
    ----------
    n : int
        approximate number of requested surface points

    Returns
    -------
    ndarray
        numpy array of xyz coordinates of surface points
    """

    u = []
    eps = 1e-10
    nequat = int(np.sqrt(np.pi*n))
    nvert = int(nequat/2)
    nu = 0
    for i in range(nvert+1):
        fi = np.pi*i/nvert
        z = np.cos(fi)
        xy = np.sin(fi)
        nhor = int(nequat*xy+eps)
        if nhor < 1:
            nhor = 1
        for j in range(nhor):
            fj = 2*np.pi*j/nhor
            x = np.cos(fj)*xy
            y = np.sin(fj)*xy
            if nu >= n:
                return np.array(u)
            nu += 1
            u.append([x, y, z])
    return np.array(u)

def vdw_surface_gamess(coordinates, elements, scale_factor, density, input_radii):
    """Computes points on the van der Waals surface of molecules.
    Parameters
    ----------
    coordinates : ndarray
        cartesian coordinates of the nuclei, in units of angstrom
    elements : list
        The symbols (e.g. C, H) for the atoms
    scale_factor : float
        The points on the molecular surface are set at a distance of
        scale_factor * vdw_radius away from each of the atoms.
    density : float
        The (approximate) number of points to generate per square angstrom
        of surface area. 1.0 is the default recommended by Kollman & Singh.
    input_radii : dict
        dictionary of user's defined VDW radii
    Returns
    -------
    radii : dict
        A dictionary of scaled VDW radii
    surface_points : ndarray
        array of the coordinates of the points on the surface
    """
    radii = {}
    surface_points = []
    # scale radii
    for i in elements:
        if i in radii.keys():
            continue
        if i in input_radii.keys():
            radii[i] = input_radii[i] * scale_factor
        elif i in vdw_r.keys():
            radii[i] = vdw_r[i] * scale_factor
        else:
            raise KeyError('%s is not a supported element; ' %i
                         + 'use the "VDW_RADII" option to add '
                         + 'its van der Waals radius.')

    # loop over atomic coordinates
    for i in range(len(coordinates)):

        # calculate approximate number of ESP grid points
        n_points = int(density * 4.0 * np.pi* np.power(radii[elements[i]], 2))

        # generate an array of n_points in a unit sphere around the atom
        dots = surface(n_points)
        #dots = fibonacci_sphere(n_points)

        # scale the unit sphere by the VDW radius and translate
        dots = coordinates[i] + radii[elements[i]] * dots

        for j in range(len(dots)):
            save = True
            for k in range(len(coordinates)):
                if i == k:
                    continue
                # exclude points within the scaled VDW radius of other atoms
                d = np.linalg.norm(dots[j] - coordinates[k])
                if d < radii[elements[k]]:
                    save = False
                    break
            if save:
                surface_points.append(dots[j])

    return np.array(surface_points), radii

def vdw_spheres(coordinates, elements, scale_factor, density, radii):

    surface_points = []
    # scale radii
    #for i in elements:
    #    if i in radii.keys():
    #        continue
    #    if i in input_radii.keys():
    #        radii[i] = input_radii[i] * scale_factor
    #    elif i in vdw_r.keys():
    #        radii[i] = vdw_r[i] * scale_factor
    #    else:
    #        raise KeyError('%s is not a supported element; ' %i
    #                     + 'use the "VDW_RADII" option to add '
    #                     + 'its van der Waals radius.')

    for i in range(len(radii)):
        radii[i] = radii[i] * scale_factor

    surfpoints = []
    for i in range(len(coordinates)):
        #n_points = int(density * 4.0 * np.pi* np.power(radii[elements[i]], 2))
        n_points = int(density * 4.0 * np.pi* np.power(radii[i], 2))
        dots = fibonacci_sphere(n_points)
        #dots = coordinates[i] + radii[elements[i]] * dots
        dots = coordinates[i] + radii[i] * dots
        surfpoints.append(dots)

    return np.array(surfpoints)

if __name__ == "__main__":


    pqr = '/home/leonardo/Desktop/Pipelines/ESP_ACCURATE/pro_esp/arch3/coulombic_grid/ar3_protein.fixed.fixed.pqr'
    llidx = pqr_scraper.find_pqr_linker_lys(pqr)
    coords = np.array(pqr_scraper.extract_ret_pqr(pqr,llidx)).T
    elements = ["C"]*len(coords)
    o1, r1 = coords[0], 1.50
    o2, r2 = coords[1], 1.50
    o3, r3 = coords[2], 1.50
    o4, r4 = coords[3], 1.50

    s1 = r1*fibonacci_sphere(10000) + o1
    s2 = r2*fibonacci_sphere(10000) + o2
    s3 = r3*fibonacci_sphere(10000) + o3
    s4 = r4*fibonacci_sphere(10000) + o4

    nucleii  = [o1, o2, o3, o4]
    vdwradii = [r1, r2, r3, r4]

    circles = multisphereintersec(nucleii, vdwradii)
    file_io.writerowsxyz('sphere1.xyz', s1)
    file_io.writerowsxyz('sphere2.xyz', s2)
    file_io.writerowsxyz('sphere3.xyz', s3)
    file_io.writerowsxyz('sphere4.xyz', s4)
    file_io.writerowsxyz('circle3.xyz', circles)
    sys.exit()

    pts = spheres_intersec(o1, o2, r1, r2)
    file_io.writerowsxyz('circle2.xyz', pts)

    start = time.time()
    surf = vdw_spheres(coords, elements, 1, 5000, vdw_r)[0]
    end = time.time()
    print(end - start)
