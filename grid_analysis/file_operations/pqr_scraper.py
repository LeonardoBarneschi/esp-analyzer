#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import re
import sys

from grid_analysis.geom_operations import molecular_linalg

#chromophore = 'RET'
#ref_atom = 'C15'

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
       ret_lys[key] = molecular_linalg.two_points_dist(a, b)

    ret_lys_sorted= sorted(ret_lys.items(), key=lambda x:x[1])
    lynker_lys = ret_lys_sorted[0][0]

    return lynker_lys

def extract_ret_pqr(pqr, llidx, heavy=None, surf=None):

    if not os.path.exists(pqr):
        return False

    llys = ['NZ', 'HZ1', 'CE', 'HE3', 'HE2']
    elements = []
    retcoord = []
    vdwcoord = []
    with open(pqr, 'r') as p:
        for line in p:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                ls = line.split()
                if len(ls) == 10 or len(ls) ==11:

                    if len(ls) == 11:
                        index = ls[5]
                        x, y, z = ls[6], ls[7], ls[8]

                    elif len(ls) == 10:
                        index = ls[4]
                        x, y, z = ls[5], ls[6], ls[7]
                if not heavy:
                    if ls[3] == 'RET':
                        elements.append(ls[2])
                        retcoord.append([x, y, z])

                    if ls[3] == 'LYS' and index == str(llidx) and ls[2] in llys:
                        elements.append(ls[2])
                        retcoord.append([x, y, z])
                else:
                    if ls[3] == 'RET' and ls[2][0] != 'H':
                        elements.append(ls[2])
                        retcoord.append([x, y, z])

                    if ls[3] == 'LYS' and index == str(llidx) and ls[2] in llys and (ls[2] == llys[0] or ls[2] == llys[2]):
                        elements.append(ls[2])
                        retcoord.append([x, y, z])

    x = [float(retcoord[i][0]) for i in range(len(retcoord))]
    y = [float(retcoord[k][1]) for k in range(len(retcoord))]
    z = [float(retcoord[j][2]) for j in range(len(retcoord))]

    if surf:
        return x, y, z, elements
    else:
        return x, y, z

def extract_ret_pi_pqr(pqr, llidx):

    pi = ['C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','NZ']

    retpi = []
    with open(pqr, 'r') as p:
        for line in p:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                ls = line.split()
                if len(ls) == 10 or len(ls) ==11:

                    if len(ls) == 11:
                        index = ls[5]
                        x, y, z = ls[6], ls[7], ls[8]

                    elif len(ls) == 10:
                        index = ls[4]
                        x, y, z = ls[5], ls[6], ls[7]

                if ls[3] == 'RET' and ls[2] in pi:
                    retpi.append([x, y, z])

                if ls[3] == 'LYS' and ls[4] == str(llidx) and ls[2] in pi:
                    retpi.append([x, y, z])

    x = [float(retpi[i][0]) for i in range(len(retpi))]
    y = [float(retpi[k][1]) for k in range(len(retpi))]
    z = [float(retpi[j][2]) for j in range(len(retpi))]

    return x, y, z

def chromophorepqr(pqr, heavy=None):

    if not os.path.exists(pqr):
        return False

    llidx = find_pqr_linker_lys(pqr)
    llys = ['NZ', 'HZ1', 'CE', 'HE3', 'HE2']
    retcoord = []
    with open(pqr, 'r') as p:

        for line in p:

            if line.startswith('ATOM') or line.startswith('HETATM'):
                ls = line.split()

                if len(ls) == 10 or len(ls) ==11:
                    if len(ls) == 11:
                        index = ls[5]
                        x, y, z = ls[6], ls[7], ls[8]

                    elif len(ls) == 10:
                        index = ls[4]
                        x, y, z = ls[5], ls[6], ls[7]

                if not heavy:
                    if ls[3] == 'RET':
                        retcoord.append(line)

                    if ls[3] == 'LYS' and ls[4] == str(llidx) and ls[2] in llys:
                        retcoord.append(line)
                else:
                    if ls[3] == 'RET' and ls[2][0] != 'H':
                        retcoord.append(line)

                    if ls[3] == 'LYS' and index == str(llidx) and ls[2] in llys and (ls[2] == llys[0] or ls[2] == llys[2]):
                        retcoord.append(line)

    return retcoord

if __name__ == "__main__":

    pqr = sys.argv[1]
    idx = find_pqr_linker_lys(pqr)
    ret = extract_ret_pqr(pqr, idx, vdw=True)
    print(ret)
