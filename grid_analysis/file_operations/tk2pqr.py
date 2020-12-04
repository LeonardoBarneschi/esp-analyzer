#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import shutil
import numpy as np
import MDAnalysis as mda
from itertools import groupby
from grid_analysis.file_operations.tk2pqr_map import standard
from grid_analysis.file_operations.tk2pqr_map import ctermin
from grid_analysis.file_operations.tk2pqr_map import ntermin
from grid_analysis.file_operations.renumber import main as renum
from grid_analysis.file_operations.parse_molcas import parse_caspt2
from grid_analysis.file_operations import getmap
import warnings

if not sys.warnoptions:
        warnings.simplefilter("ignore")

pqr_neutrans = {'CD':  0.000000,
                'HD3': 0.062100,
                'HD2': 0.062100,
                'CG':  0.018700,
                'HG3': 0.010300,
                'HG2': 0.010300,
                'CB': -0.009400,
                'HB3': 0.036200,
                'HB2': 0.036200,
                'CA': -0.240000,
                'HA':  0.142600,
                'N':  -0.398050,
                'H':   0.224550,
                'C':   0.683950,
                'O':  -0.639550,
                'HZ1': 0.000000,
                'NZ':  0.000000,
                'CE':  0.000000,
                'HE2': 0.000000,
                'HE3': 0.000000}

def two_points_dist(a, b):
    c = a - b
    return np.linalg.norm(c)

def find_pqr_linker_lys(pdb, ref_atom='C15', chromophore='RET'):

    lysines = {}
    ret_lys = {}
    ref = []

    if not os.path.exists(pdb):
        print('The file {0} was not found'.format(pdb))
        return False

    with open(pdb , 'r') as p:
        for line in p:
            ls = line.split()
            if len(ls) == 10 or len(ls) ==11:

                if len(ls) == 11:
                    index = ls[5]
                    x, y, z = ls[6], ls[7], ls[8]

                elif len(ls) == 10:
                    index = ls[4]
                    x, y, z = ls[5], ls[6], ls[7]

                if ls[3] == chromophore and ls[2] == ref_atom:
                    ref.append([x,y,z])
                    ref = ref[0] # flatten list

                if ls[3] == 'LYS' and ls[2] == 'NZ':
                    lysines[index] = np.asarray([float(x), float(y), float(z)])
    # DEBUG
    #print("\n".join("{}\t{}".format(k, v) for k, v in lysines.items()))

    ref = np.asarray([float(k) for k in ref])
    for key in lysines.keys():
       a = ref
       b = lysines[key]
       ret_lys[key] = two_points_dist(a, b)

    ret_lys_sorted= sorted(ret_lys.items(), key=lambda x:x[1])
    lynker_lys = ret_lys_sorted[0][0]

    return lynker_lys

def fix_pqr_linkerlys(pqr, lys, out):

    if not os.path.exists(pqr):
        print('The file {0} was not found'.format(pqr))
        return False

    pqr_fixed = []
    with open(pqr,'r') as p:
        for line in p:
            ls = line.split()

            if ls[5] == lys and ls[2] in pqr_neutrans.keys():
                line = line[:58]+'%7.4f' % pqr_neutrans[ls[2]]+line[65:]

            pqr_fixed.append(line)

    #if not os.path.exists(out):

    with open(out, 'w+') as o:
        for l in pqr_fixed:
            o.write(l)

    return

def fixpdb(oldpdb, newpdb):
    newpdblist = []
    if os.path.exists(oldpdb):
        with open(oldpdb, 'r') as oldpdb:
            for line in oldpdb:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    line = line[:21] + 'A' + line[22:]
                    newpdblist.append(line)

    #if not os.path.exists(newpdb):
    with open(newpdb, 'w') as np:
        np.writelines(newpdblist)

    return newpdb

def tk2pdb(oldpdb, oldxyz, newxyz, outpdb="new.pdb"):
    '''
    Function to convert a Tinker XYZ file to PDB. The starting Tinker XYZ and
    PDB file used to generate it are needed to obtain a map for the conversion
    of atom indices.

    Parameters
    ----------
    oldpdb: str.
        Initial PDB file name.
    oldxyz: str.
        Initial Tinker XYZ file name.
    newxyz: str.
        Tinker XYZ file to be converted to PDB.
    outpdb: str.
        Output PDB file name.
    '''

    # Get old pdb coordinates
    oldp = mda.Universe(oldpdb)
    oldpcoords = oldp.select_atoms("all").positions

    # Get old xyz coordinates, except link atom
    oldx = mda.Universe(oldxyz, topology_format="TXYZ")
    oldxcoords = oldx.select_atoms("all").positions

    # Get new xyz coordinates, except link atom
    newx = mda.Universe(newxyz, topology_format="TXYZ")
    newxcoords = newx.select_atoms("all").positions[:-1]

    # Get the map from old to new
    # We do coordinate comparison in fortran for performance 
    old_to_new = getmap.getmap(oldpcoords, oldxcoords)

    # Convert fortran idxs to python idxs
    oldidxs = old_to_new[:,0].astype(int) - 1
    newidxs = old_to_new[:,1].astype(int) - 1

    # Reorder coordinates and print
    newp = mda.Universe(oldpdb)
    newp.select_atoms("all").positions = newxcoords[newidxs]
    newp.select_atoms("all").write(outpdb)

    return outpdb

def tk2pqr(oldpdb, oldxyz, newxyz, output=False, chrm="RET"):
    '''
    Convert a Tinker-xyz file to PQR format.
    - Phase 1) Create a map (indices) using corresponding
               TKXYZ and PDB files.
    - Phase 2)
    - Phase 3)

    Parameters
    ----------
    oldpdb : file, str, or os.PathLike
           input PDB. Required for indices mapping.

    oldxyz : file, str, or os.PathLike
           input TKXYZ. Required for indices mapping.

    newpdb : file, str, or os.PathLike

    output : file, str, or os.PathLike

    Returns
    -------
    '''

    if not output:
        newpdb = oldpdb
        newpqr = newxyz[:-4]+'.pqr'
        pqrs0 =  newxyz[:-4]+'_s0.pqr'
        pqrs1 =  newxyz[:-4]+'_s1.pqr'

    else:
        newpdb = newxyz[:-4]+'_new.pdb'
        newpqr = newxyz[:-4]+'_new.pqr'
        pqrs0 =  output[:-4]+'_s0.pqr'
        pqrs1 =  output[:-4]+'_s1.pqr'

    a = mda.Universe(oldpdb)
    starting_resid = a.atoms.resids[0]
    renum(oldpdb, newpdb, starting_resid)
    newpdb = fixpdb(newpdb, newpdb)
    outpdb = tk2pdb(newpdb, oldxyz, newxyz, newpqr)

    u = mda.Universe(newpqr)

    start = None

    # TODO -> fix for GAPS in sequence
    np.where(u.atoms.resnames == chrm)
    ref = np.where(u.atoms.resnames == chrm)[0][0]
    end = u.atoms.resids[ref-1]

    charges = []
    vdwradii = []

    for i in range(u.atoms.names.shape[0]):

        resname = u.atoms.resnames[i]
        atname = u.atoms.names[i]
        resid = u.atoms.resids[i]
        mapping_dict = standard[resname].copy()

        if not start or start == resid:
            mapping_dict.update(ntermin[resname])
            start = resid

        if resid == end:
            mapping_dict.update(ctermin[resname])

        try:
            charge = mapping_dict[atname]["charge"]
            vdwrad = mapping_dict[atname]["vdw"]
            charges.append(charge)
            vdwradii.append(vdwrad)
        except:
            print('problem with: {0} | {1} | {2}\n'.format(resname, atname, resid))

    u.atoms.charges = np.array(charges)
    u.atoms.radii = np.array(vdwradii)


    a = np.c_[u.atoms.resids, np.arange(u.atoms.resids.shape[0])]
    groups = np.split(a[:,1], np.unique(a[:, 0], return_index=True)[1][1:])

    for group in groups:
        resid = u.atoms.resids[group[0]]
        resn = u.atoms.resnames[group[0]]

    # Create Opsin PQR file (keep both writes)
    u.select_atoms("all").write(newpqr)
    lab = find_pqr_linker_lys(newpqr)
    fix_pqr_linkerlys(newpqr, lab, newpqr)
    u = mda.Universe(newpqr)
    print('\n Charge check for input: {} --> {:.3f}'.format(newxyz, u.atoms.charges.sum()))
    u.select_atoms("all").write(newpqr)

    # Create Chromophore PQR files
    # TODO GENERALIZE (ESPECIALLY ESPF CHARGES)
    atoms, structure, _, qs = parse_caspt2(output)
    atoms = atoms[-54:]
    qs = qs[-54:]
    structure = structure[-54:]
    pqrret = getmap.getmap(structure, u.atoms.positions)-1
    pqrret = pqrret.astype(int)[:-1]
    print(pqrret)
    #sys.exit()
    idxs = u.select_atoms("name CE HE2 HE3 and resid %s" % lab).indices
    retidxs = []
    for i, idx in enumerate(pqrret[:,1]):
        if idx in idxs:
            retidxs.append(i)
    retidxs = np.array(retidxs)
    qs[retidxs] += qs[-1]/len(retidxs)

    # do not take LAH into account
    qs = qs[:-1]

    # set opsin charges to zero
    u.atoms.charges = np.zeros_like(u.atoms.charges)
    tmpqs = np.zeros_like(u.atoms.charges)

    # GS (but can be generalized easily)
    tmpqs[pqrret[:,1]] = qs[pqrret[:,0],0]
    u.atoms.charges = tmpqs
    print('\n Charge check for state: {} --> {:.3f}'.format('S0', u.atoms.charges.sum()))
    u.select_atoms("all").write(pqrs0)

    # ES
    tmpqs[pqrret[:,1]] = qs[pqrret[:,0],1]
    u.atoms.charges = tmpqs
    print('\n Charge check for state: {} --> {:.3f}'.format('S1', u.atoms.charges.sum()))
    u.select_atoms("all").write(pqrs1)

    #return newpqr, pqrs0, pqrs1
    return newpqr, pqrs1

if __name__ == "__main__":

    pass

