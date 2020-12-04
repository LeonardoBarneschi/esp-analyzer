#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import subprocess
import sys

def fix_pqr_mwfn(pqr, tmppqr=None):
    """
    Multiwfn required format (python idxs):
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    """

    if not tmppqr:
        tmppqr = pqr[:-3]+'tmp.pqr'

    fmt = "%-6s%5s %4s %3s %1s %3s     %7s %7s %7s %8s  %6s\n"
    newlines = []
    with open(pqr, 'r') as p:
        for line in p:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                newline = fmt % tuple(line.split())
                newlines.append(newline)

    with open(tmppqr, 'w') as p:
        for l in newlines:
            p.write(l)

    return tmppqr

def calculate_esp(pqr, wfngrid, out=None):
    '''
    Coulombic ESP calculation by Multiwfn

    Parameters
    ----------
    pqr : file, str, or os.PathLike
        PQR file (standard format)

    wfn : file, str, or os.PathLike
        Communicator file for Multiwfn. Specifies the grid
        points for ESP calculation. The first row specifies
        the n. of grid points, the rest of the file is
        formatted as x, y, z coordinates.

    Returns
    -------
    df : a Pandas dataframe with 4-columns
       (x, y, z, esp)
    '''

    pqr = fix_pqr_mwfn(pqr)
    args = locals()
    argout = args['out']

    if not out:
        out = 'communicator'

    if not (sys.version_info.major == 3 and sys.version_info.minor >= 5):
        cmd = "echo '{0}\n 5\n 8\n 100\n {1}\n {2}\n' | Multiwfn".format(pqr, wfngrid, out)
        result = subprocess.call([cmd], shell=True, capture_output=True)

    else:
        cmd = "echo '{0}\n 5\n 8\n 100\n {1}\n {2}\n' | Multiwfn".format(pqr, wfngrid, out)
        result = subprocess.run([cmd], shell=True, capture_output=True)

    return out

def esp_from_edens(orb, grid, out=None):
    if not out:
        out = 'communicator'

    if not (sys.version_info.major == 3 and sys.version_info.minor >= 5):
        cmd = "echo '{0}\n\n 5\n 12\n 100\n {1}\n {2}\n' | Multiwfn".format(orb, grid, out)
        result = subprocess.call([cmd], shell=True, capture_output=True)

    else:
        cmd = "echo '{0}\n\n 5\n 12\n 100\n {1}\n {2}\n' | Multiwfn".format(orb, grid, out)
        result = subprocess.run([cmd], shell=True, capture_output=True)

    return out

def edens2cub(orb):

    if not (sys.version_info.major == 3 and sys.version_info.minor >= 5):
        cmd = "echo '{0}\n\n 5\n 1\n 1\n 2\n' | Multiwfn".format(orb)
        result = subprocess.call([cmd], shell=True, capture_output=True)

    else:
        cmd = "echo '{0}\n\n 5\n 1\n 1\n 2\n' | Multiwfn".format(orb)
        result = subprocess.run([cmd], shell=True, capture_output=True)

    return 'density.cub'

if __name__ == "__main__":

    pqr = sys.argv[1]
    grd = sys.argv[2]
    out = sys.argv[3]
    calculate_esp(pqr, grd, out)
