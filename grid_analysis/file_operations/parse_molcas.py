#!/usr/bin/env python

import os
import re
import sys
import csv
import numpy as np
import pandas as pd
import MDAnalysis as mda
import multiprocessing as mp

au2eV = 27.21138505
au2fs = 41.341374575751
au2kcalmol = 627.50


def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)


def parse_caspt2(filename):
    '''
    Function to parse a Molcas CASPT2 calculation output.

    Parameters
    ----------
    filename: str.
        Name of the file to parse.

    Returns
    -------
    atoms: list.
        list of atom symbols.
    structure: np.array (N, 3).
        array of atomic coordinates (in angstroem).
    energies: np.array (M).
        array of CASPT2 energies. The n-th index corresponds to the n-th state.
    foscs: np.array (M, M).
        array of the intensity of transition from state A to state B.
    occs: np.array (M, P).
        array of Natural Orbital occupancies for the orbitals in the active
        space.
    '''

    with open(filename) as f:

        atoms = []
        structure = []
        energies = []
        foscs = []
        idxs = []
        occs = []
        chgs = []
        espf = []
        nstates = []

        for line in f:

            if "Cartesian coordinates in Angstrom:".lower() in line.lower():
                line = skiplines(f, 3)
                while len(line.split()) == 5:
                    data = line.split()
                    atom = data[1]
                    coords = list(map(float, data[-3:])) # coordinates in angstroem
                    # coords = list(map(float, data[2:5])) # coordinates in au
                    atoms.append(atom)
                    structure.append(coords)
                    line = next(f)

            if "N-E" in line:
                data = list(map(float, line.split()[1:]))
                chgs.extend(data)

            if "Mulliken population analysis for root number:".lower() in line.lower():
                nstate = int(line.split()[-1])
                nstates.append(nstate)

            if "Charge on" in line:
                q = float(line.split()[-1])
                espf.append(q)

            if "CASPT2 Root" in line:
                data = line.split()
                ene = float(data[-1])
                energies.append(ene)

            if "Dipole transition strengths" in line:
                line = skiplines(f, 5)
                data = line.split()
                while len(data) == 7:
                    idx1 = int(data[0])
                    idx2 = int(data[1])
                    fosc = float(data[2])
                    idxs.append(tuple([ idx1, idx2 ]))
                    foscs.append(fosc)
                    line = skiplines(f)
                    data = line.split()

            if "Natural orbitals and occupation numbers" in line:
                line = skiplines(f)
                stateoccs = []
                while len(line.strip()) > 0:
                    data = line.split()
                    try:
                        tmpoccs = list(map(float, data))
                    except ValueError:
                        tmpoccs = list(map(float, data[2:]))

                    stateoccs.extend(tmpoccs)
                    line = skiplines(f)

                occs.append(stateoccs)

        try:
            nstates = max(nstates)
            atoms = np.array(atoms)
            structure = np.array(structure)
            #energies = np.array(energies)
            #energies -= energies.min()
            #energies *= au2eV
            #occs = np.array(occs)
            chgs = np.array(chgs).reshape(nstates, -1).T
            espf = np.array(espf).reshape(nstates, -1).T
            # Create a matrix with all foscs printed from Molcas
            #M = np.zeros( (len(energies), len(energies)) )
            #for n, idx in enumerate(idxs):
            #    i = idx[0] - 1
            #    j = idx[1] - 1
            #    M[i,j] = foscs[n]

            #return atoms, structure, energies, M, occs, chgs
            return atoms, structure, chgs, espf

        except:
            return None, None, None, None


def test_occs(occs, hb=1.98, lb=0.02):
    '''
    Function to test Natural Orbital occupancies. If the occupancy for all
    orbitals in the active space is between the bounds lb and hb, then this
    returns False, and the active space is fine.

    Parameters
    ----------
    occs: np.array (M, P).
        array of Natural Orbital occupancies for the orbitals in the active
        space.
    hb: float (default: 1.98).
        higher bound for NBO occupancy.
    lb: float (default: 0.02).
        lower bound for NBO occupancy.

    Returns
    -------
    warning: bool.
        whether the active space for the input occupancies is to be checked.
    '''

    test1 = occs[occs >= hb]
    test2 = occs[occs <= lb]

    if np.any([ test1, test2 ]):
        warning = True
    else:
        warning = False

    return warning


def ene_trj(filename, nstates=2):

    cvt = lambda x: float(x.replace("D", "e"))
    cvt_d = { i : cvt for i in range(4 + nstates) }
    data = np.loadtxt(filename, skiprows=1, converters=cvt_d, encoding=None)
    times = data[:,0] / au2fs
    enes = data[:,-nstates:]

    return times, enes


def geom_trj(top, traj, ICs):

    trj = mda.Universe(top, traj)
    nframes = trj.trajectory.n_frames
    nICs = len(ICs)
    data = np.zeros((nframes, nICs))
    for i, ts in enumerate(trj.trajectory):
        for j, IC in enumerate(ICs):

            nats = len(IC)
            if nats == 2:
                A, B = ts.positions[ [ IC ]]
                v = np.linalg.norm(A - B)

            if nats == 3:
                A, B, C = ts.positions[ [ IC ]]
                v = angle(A, B, C)

            if nats == 4:
                A, B, C, D = ts.positions[ [ IC ]]
                v = dihedral(A, B, C, D)

            data[i,j] = v

    return data


def analyse_trj(fname, ICs, nstates):

    enename = fname[:-3] + "md.energies"
    trjname = fname[:-3] + "md.xyz"

    times, enes = ene_trj(enename, nstates=nstates)

    data = geom_trj(trjname, trjname, ICs)
    res = np.c_[ enes, data ]
    df = pd.DataFrame(index=times, data=res)

    return df


def parallel_fn(fn, iterable, nproc=None):
    '''
    Function to execute a generic function in parallel applying it to all
    the elements of an iterable.

    Parameters
    ----------
    fn: object.
        Python standalone function.
    iterable: iterable.
        Collection of elements on which fn should be applied.
    nproc: integer (default: None).
        Number of processors for the parallelisation. If not specified all
        available processors will be used.

    Returns
    -------
    data: list.
        List of the returns of function fn for each element of iterable.
    '''

    if not nproc:
        nproc = os.cpu_count()

    pool = mp.Pool(nproc)
    data = pool.map(fn, iterable)
    pool.close()
    pool.join()

    return data


def v1v2_angle(v1, v2):
    '''
    Function to compute the angle between two vectors.

    Parameters
    ----------
    v1: np.array (3).
        First vector.
    v2: np.array (3).
        Second vector.

    Returns
    -------
    theta: float.
        Angle between the vector (in degrees).
    '''

    try:
        theta = np.degrees(np.arccos(
                    np.dot(v1, v2) / ( np.linalg.norm(v1) * np.linalg.norm(v2) )
                ))
    except:
        theta = 0.0

    return theta


def angle(A, B, C):
    '''
    Function to compute the angle defined by points A, B, and C.

    Parameters
    ----------
    A: np.array (2 or 3).
        First point.
    B: np.array (2 or 3).
        Second point, vertex of the angle.
    C: np.array (2 or 3).
        Third point.

    Returns
    -------
    abc: float.
        Angle ABC (in degrees).
    '''

    ab = B - A
    bc = C - B
    abc = v1v2_angle(ab, bc)

    return abc


def dihedral(A, B, C, D):
    '''
    Function to compute the dihedral angle between the planes containing
    segments AB and CD.

    Parameters
    ----------
    A: np.array (3).
        First point.
    B: np.array (3).
        Second point.
    C: np.array (3).
        Third point.
    D: np.array (3).
        Fourth point.

    Returns
    -------
    dihedral: float.
        Dihedral angle defined by AB and CD (in degrees).
    '''

    b0 = -(B - A)
    b1 = C - B
    b2 = D - C

    # normalize b1 so that it does not influence magnitude of vector
    # projections that come next
    b1 /= np.linalg.norm(b1)

    # vector projections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    dihedral = np.arctan2(y, x) * 180 / np.pi

    return dihedral


if __name__ == '__main__':
    pass
