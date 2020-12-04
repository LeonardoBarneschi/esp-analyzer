#!/usr/bin/env python3

import sys
import numpy as np

ang2au = 1/(0.52917721067121)
au2ang = 0.52917721067121

def grid_ang2au(grid, ovw=False):

    out = grid[:-4]+'_new.dat'

    xyz = np.loadtxt(grid, skiprows=1)
    xyz = xyz * ang2au

    if not ovw:
        with open(grid[:-4]+'_new.dat', 'w') as g:
            g.write(str(xyz.shape[0])+'\n')
            np.savetxt(g, xyz, fmt="%10.6f")

        return out

    else:
        print(grid)
        with open(grid, 'w') as go:
            go.write(str(xyz.shape[0])+'\n')
            np.savetxt(go, xyz, fmt="%10.6f")

        return grid

def grid_au2ang(grid, ovw=False):

    out = grid[:-4]+'_new.dat'

    xyz = np.loadtxt(grid, skiprows=1)
    xyz = xyz * ang2au

    if not ovw:
        with open(grid[:-4]+'_new.dat', 'w') as g:
            g.write(str(xyz.shape[0])+'\n')
            np.savetxt(g, xyz, fmt="%10.6f")

        return out

    else:
        print(grid)
        with open(grid, 'w') as go:
            go.write(str(xyz.shape[0])+'\n')
            np.savetxt(go, xyz, fmt="%10.6f")

        return grid

if __name__ == "__main__":
    grid = sys.argv[1]
    grid_ang2au(grid)
