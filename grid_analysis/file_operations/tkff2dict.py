#!/usr/bin/env python3

import re
import sys
import glob
import json
import shutil
import numpy as np
import MDAnalysis as mda
from collections import OrderedDict

NAA = {"N-Term ALA"            : "ALA",
       "N-Term ARG"            : "ARG",
       "N-Term ASN"            : "ASN",
       "N-Term ASP"            : "ASP",
       "N-Term CYS (-SH)"      : "CYS",
       "N-Term CYS (-SS-)"     : "CYS",
       "N-Term GLN"            : "GLN",
       "N-Term GLU"            : "GLU",
       "N-Term GLY"            : "GLY",
       "N-Term HIS (+)"        : "HIS",
       "N-Term HIS (HD)"       : "HID",
       "N-Term HIS (HE)"       : "HIE",
       "N-Term ILE"            : "ILE",
       "N-Term LEU"            : "LEU",
       "N-Term LYS"            : "LYS",
       "N-Term MET"            : "MET",
       "N-Term PHE"            : "PHE",
       "N-Term PRO"            : "PRO",
       "N-Term SER"            : "SER",
       "N-Term THR"            : "THR",
       "N-Term TRP"            : "TRP",
       "N-Term TYR"            : "TYR",
       "N-Term VAL"            : "VAL"}

CAA = {"C-Term ALA"            : "ALA",
       "C-Term ARG"            : "ARG",
       "C-Term ASN"            : "ASN",
       "C-Term ASP"            : "ASP",
       "C-Term CYS (-SH)"      : "CYS",
       "C-Term CYX (-SS-)"     : "CYS",
       "C-Term GLN"            : "GLN",
       "C-Term GLU"            : "GLU",
       "C-Term GLY"            : "GLY",
       "C-Term HIS (+)"        : "HIS",
       "C-Term HIS (HD)"       : "HID",
       "C-Term HIS (HE)"       : "HIE",
       "C-Term ILE"            : "ILE",
       "C-Term LEU"            : "LEU",
       "C-Term LYS"            : "LYS",
       "C-Term MET"            : "MET",
       "C-Term PHE"            : "PHE",
       "C-Term PRO"            : "PRO",
       "C-Term SER"            : "SER",
       "C-Term THR"            : "THR",
       "C-Term TRP"            : "TRP",
       "C-Term TYR"            : "TYR",
       "C-Term VAL"            : "VAL"}

SAA = {"ALANINE"               : "ALA",
       "ARGININE"              : "ARG",
       "ASPARAGINE"            : "ASN",
       "ASPARTIC ACID"         : "ASP",
       "ASPARTIC ACID NEUTRAL" : "ASH",
       "CYSTEINE (-SH)"        : "CYS",
       "CYSTINE (-SS-)"        : "CYS",
       "GLUTAMINE"             : "GLN",
       "GLUTAMIC ACID"         : "GLU",
       "GLUTAMIC ACID NEUTRAL" : "GLH",
       "GLYCINE"               : "GLY",
       "HISTIDINE (+)"         : "HIS",
       "HISTIDINE (HD)"        : "HID",
       "HISTIDINE (HE)"        : "HIE",
       "ISOLEUCINE"            : "ILE",
       "LEUCINE"               : "LEU",
       "LYSINE"                : "LYS",
       "LYSINE NEUTRAL"        : "LYD",
       "METHIONINE"            : "MET",
       "PHENYLALANINE"         : "PHE",
       "PROLINE"               : "PRO",
       "SERINE"                : "SER",
       "THREONINE"             : "THR",
       "TRYPTOPHAN"            : "TRP",
       "TYROSINE "             : "TYR",
       "VALINE"                : "VAL"}

# You have to wish strongly that CX-CN actually 
# correspond to sp2 (C2R) and sp3 (C3R) carbons.
# Depends really on how the PQR is written
NSAA = {"C2R" : "RET",
        "C3R" : "RET",
        "HR"  : "RET",
        "CL-" : "CL",
        "NA+" : "NA",
        "OT"  : "OT",
        "HT"  : "HT"}

TKXYZPQR = {'C2R' : ('C5','C6','C7','C8','C9','C10',
                     'C11','C12','C13','C14','C15'),
            'C3R' : ('C1','C2','C3','C4','C16','C17','C18','C19','C20'),
            'HR'  : ('H21','H22','H23','H24','H25','H26','H27',
                     'H28','H29','H30','H31','H32','H33','H33',
                     'H34','H35','H36','H37','H38','H39','H40',
                     'H41','H42','H43','H44','H45','H46','H47',
                     'H48'),
            'CL-' : ('CL'),
            'NA+' : ('NA'),
            'OT' : ('OT'),
            'HT' : ('HT')}

def name2resname(ffdict):
    '''
    The Tinker formatter FF file (AMBER94) + Chromophores
    here called melacu63.prm is read line by line in a dictionary
    where values are the whole line strings.
    '''
    ctermin = {}
    ntermin = {}
    standard = {}
    nonstandard = {}

    for i, line in enumerate(ffdict["atom"], start=1):
        atomname = line.split('"')[1].split()[-1]
        atomtype  = "type "+line.split()[1]
        atomclass = "class "+line.split()[2]
        for k in SAA.keys():
            if '"'+k.lower() in line.lower():
                try:
                   standard[SAA[k]][atomname] = {}
                   standard[SAA[k]][atomname][atomtype] = {}
                   standard[SAA[k]][atomname][atomclass] = {}
                except KeyError:
                   standard[SAA[k]] = {}
                   standard[SAA[k]][atomname] = {}
                   standard[SAA[k]][atomname][atomtype] = {}
                   standard[SAA[k]][atomname][atomclass] = {}

        for w in NAA.keys():
            if w.lower() in line.lower():
                try:
                    ntermin[NAA[w]][atomname] = {}
                    ntermin[NAA[w]][atomname][atomtype] = {}
                    ntermin[NAA[w]][atomname][atomclass] = {}
                except KeyError:
                    ntermin[NAA[w]] = {}
                    ntermin[NAA[w]][atomname] = {}
                    ntermin[NAA[w]][atomname][atomtype] = {}
                    ntermin[NAA[w]][atomname][atomclass] = {}
        for u in CAA.keys():
            if u.lower() in line.lower():
                try:
                    ctermin[CAA[u]][atomname] = {}
                    ctermin[CAA[u]][atomname][atomtype] = {}
                    ctermin[CAA[u]][atomname][atomclass] = {}
                except KeyError:
                    ctermin[CAA[u]] = {}
                    ctermin[CAA[u]][atomname] = {}
                    ctermin[CAA[u]][atomname][atomtype] = {}
                    ctermin[CAA[u]][atomname][atomclass] = {}

        for j in NSAA.keys():
            if j.upper() in line.split()[3].upper():
                for i in TKXYZPQR.keys():
                    if i == j:
                        try:
                            nonstandard[NSAA[j]][i] = {}
                            nonstandard[NSAA[j]][i][atomtype] = {}
                            nonstandard[NSAA[j]][i][atomclass] = {}
                        except KeyError:
                            nonstandard[NSAA[j]] = {}
                            nonstandard[NSAA[j]][i] = {}
                            nonstandard[NSAA[j]][i][atomtype] = {}
                            nonstandard[NSAA[j]][i][atomclass] = {}

    return standard, ctermin, ntermin, nonstandard

def name2vdw(ffdict, mapdict, nonstd=False):

    for line in ffdict["vdw"]:
        data = line.split()
        oldclass = data[1]
        vdwradii = float(data[2])
        for reskey in mapdict.keys():
            for namekey in mapdict[reskey].keys():
                for atomkey in mapdict[reskey][namekey].keys():
                    if 'class' in atomkey:
                        if atomkey.split()[1] == oldclass:
                            mapdict[reskey][namekey][atomkey] = vdwradii

    for line in ffdict["charge"]:
        oldattype = line.split()[1]
        charge = float(line.split()[2])
        for reskey in mapdict.keys():
            for namekey in mapdict[reskey].keys():
                for atomkey in mapdict[reskey][namekey].keys():
                    if 'type' in atomkey:
                        if atomkey.split()[1] == oldattype:
                            mapdict[reskey][namekey][atomkey] = charge

    if nonstd:
        for k in mapdict.keys():
            print(k)
            for w in mapdict[k].keys():
                for l in mapdict[k][w].keys():
                    if 'type' in l:
                        mapdict[k][w][l] = 0.000
    return mapdict

def read_tk(prmfile):

    secs = [ "forcefield",
             "vdwtype",
             "radiusrule",
             "radiustype",
             "radiussize",
             "epsilonrule",
             "vdw-14-scale",
             "chg-14-scale",
             "electric",
             "dielectric",
             "atom",
             "vdw",
             "bond",
             "angle",
             "ureybrad",
             "imptors",
             "torsion",
             "charge",
             "biotype" ]

    sections = OrderedDict()

    with open(prmfile) as f:
        for line in f:

            line = line.rstrip()
            if line.startswith("#"):
                pass

            if line.strip():
                data = line.split()
                if data[0] in secs:
                    try:
                        sections[data[0]].append(line)
                    except:
                        sections[data[0]] = [ line ]

    return sections

if __name__ == '__main__':

    ffpath = sys.argv[1]
    sec = read_tk(ffpath)
    standard, ctermin, ntermin, nonstandard = name2resname(sec)

    vdwchg1 = name2vdw(sec, standard)
    vdwchg2 = name2vdw(sec, ctermin)
    vdwchg3 = name2vdw(sec, ntermin)
    vdwchg4 = name2vdw(sec, nonstandard, nonstd=True)

    print('\n')
    print("\n".join("{}\t{}".format(k, v) for k, v in vdwchg1.items()))
    print('\n')
    print("\n".join("{}\t{}".format(k, v) for k, v in vdwchg2.items()))
    print('\n')
    print("\n".join("{}\t{}".format(k, v) for k, v in vdwchg3.items()))
    print('\n')
    print("\n".join("{}\t{}".format(k, v) for k, v in vdwchg4.items()))
    print('\n')

    with open('standard.json', 'w') as fp1:
        json.dump(vdwchg1, fp1, indent=4)
    with open('ctermin.json', 'w') as fp2:
        json.dump(vdwchg2, fp2,sort_keys=True, indent=4)
    with open('ntermin.json', 'w') as fp3:
        json.dump(vdwchg3, fp3,sort_keys=True, indent=4)
    with open('nonstandard.json', 'w') as fp4:
        json.dump(vdwchg4, fp4,sort_keys=True, indent=4)

