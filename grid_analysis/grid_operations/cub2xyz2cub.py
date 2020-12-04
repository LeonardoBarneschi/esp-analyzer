#! /usr/bin/env python3

import csv, os, sys, argparse
import pandas as pd
import numpy as np

from grid_analysis.externals import multiwfn as mwfn
from grid_analysis.geom_operations import molecular_linalg as ml
from grid_analysis.grid_operations.cubmanager import Cube
from grid_analysis.file_operations import tk2pqr

au2ang = 0.5291772083
ang2au = 1/au2ang
au2volts = 27.211386245988

def read_cube(cub):

    lidx = 0
    svec = []
    vals = []
    data = []
    geom = []
    # Read cube file and parse all data
    for line in open(cub, 'r'):
        lidx += 1
        if lidx == 3:
           try:
               nat=int(line.split()[0])
               origin=[float(line.split()[1]),
                       float(line.split()[2]),
                       float(line.split()[3])]
           except:
               print("ERROR: cube format non recognized")

        elif lidx > 3 and lidx <= 6:
            svec.append(line.split())
        elif lidx > 6 and lidx <= 6+nat:
            geom.append(line.split())
        elif lidx > 5:
            if lidx > 6+nat:
                for i in line.split():
                    vals.append(float(i))

    idx=-1

    for i in range(0,int(svec[0][0])):
        for j in range(0,int(svec[1][0])):
            for k in range(0,int(svec[2][0])):
                idx+=1
                x = origin[0] + i * float(svec[0][1])
                y = origin[1] + j * float(svec[1][2])
                z = origin[2] + k * float(svec[2][3])
                data.append([x, y, z, vals[idx]])

    nx = len(range(0,int(svec[0][0])))
    ny = len(range(0,int(svec[1][0])))
    nz = len(range(0,int(svec[2][0])))

    data = np.array(data)

    return cub, data, (nx, ny, nz)

def edens2esp(cub, data, orb, vox, iso=0.001, gridesp=None, wfn=None, calcesp=False, out=None):

    nx, ny, nz = vox[0], vox[1], vox[2]

    if not wfn:
        wfn = cub[:-4]+'.wfn'

    if not out:
        out = cub[:-4]+'_new.cub'

    rho = data[:,-1]
    #rho = np.round(data[:,-1], decimals=3)
    idxs = (rho >= 0.00092) & (rho <= 0.0018)
    #idxs = rho == 0.001
    grid = data[idxs,:-1]

    # Generate Isosurface in grid (x,y,z,rho) format
    with open(wfn,'w') as o:
        o.write("%d \n" % (len(grid)))
        np.savetxt(o, grid, fmt='%10.6f')

    if calcesp:
        # Evaluate ESP at isosurface
        gridesp = mwfn.esp_from_edens(orb, wfn)
        gridesp = np.loadtxt(gridesp, skiprows=1)
        v = np.zeros(len(data))
        v[idxs] = gridesp[:,-1]
        #v = v.reshape(nx, ny, nz)

        cub = Cube(cub)
        cub.data = v
        cub.dump(out)

    return out

def espcub(cub, data, pqr, mode=0, gridesp=None, wfn=None, surfile=None, out=None):

    if not wfn:
        wfn = cub[:-4]+'.wfn'

    if not out:
        out = cub[:-4]+'_new.cub'

    rho = data[:,-1]
    # Just take all the points
    idxs = (rho >= -10000) & (rho <= 10000)
    grid = data[idxs,:-1]

    # Generate Isosurface in grid (x,y,z,rho) format
    with open(wfn,'w') as o:
        o.write("%d \n" % (len(grid)))
        np.savetxt(o, grid, fmt='%10.6f')

    # Evaluate ESP over the whole electron density
    # CONVERT TO ANGSTROM (FROM A.U OF MULTIWFN)
    V, coords = ml.coulomb_esp(pqr, mode, surfile=wfn)

    gridesp = ml.write_coulomb_esp('tmp.esp', V, coords * ang2au, ovw=True)
    gridesp = np.loadtxt(gridesp, skiprows=1)
    os.remove('tmp.esp')
    v = np.zeros(len(data))
    v[idxs] = gridesp[:,-1]
    cub = Cube(cub)

    # Convert to Volts
    cub.data = v * au2volts
    cub.dump(out)

    return out

if __name__ == "__main__":

    #oldpdb = sys.argv[1]
    #oldxyz = sys.argv[2]
    #newxyz = sys.argv[3]
    #output = sys.argv[4]
    #pqr_pro, pqr_ret = tk2pqr(oldpdb, oldxyz, newxyz, output)
    #sys.exit()

    # From PQR
    cub = sys.argv[1]
    orb = sys.argv[2]
    pqr = sys.argv[3]
    cub, data, (nx, ny, nz) = read_cube(cub)
    cubesp = espcub(cub, data, pqr, mode=0)

    ### cub = sys.argv[1]
    ### orb = sys.argv[2]
    ### esp = sys.argv[3]
    ### cub, data, (nx, ny, nz) = read_cube(cub)
    ### cubesp = edens2esp(cub, data, orb, (nx, ny, nz), calcesp=False)#, gridesp='communicator')
