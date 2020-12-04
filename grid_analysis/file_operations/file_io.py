#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os

def writerowsxyz(filename, rows):
    '''
    '''
    print(rows.shape)
    if not os.path.exists(filename):
        try:
            with open(filename, 'w') as f:
                print('rows len is: ', len(rows))
                f.write('{0} \n\n'.format(len(rows)))
                np.savetxt(f, np.c_[np.ones(len(rows)), rows], fmt='%3d %10.6f %10.6f %10.6f')
                return filename
        except:
            print('An error occured while trying to open {0}. \n'.format(filename))
            return False

    else:
        print('\n The file {0}, already exists. \n'.format(filename))
        return False

def writerowswfn(filename, rows):
    '''
    '''
    if not os.path.exists(filename):
        try:
            with open(filename, 'w') as f:
                f.write('{0}\n'.format(len(rows)))
                np.savetxt(f, rows, fmt='%10.6f')
            return filename
        except:
            print('An error occured while trying to open {0}. \n'.format(filename))
            return False

    else:
        print('\n The file {0}, already exists. \n'.format(filename))
        return False
