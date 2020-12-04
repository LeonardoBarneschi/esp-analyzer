#!/usr/bin/env python

import numpy as np
from getpot import getpot

from scipy.spatial.distance import cdist

if __name__ == '__main__':

    ang2au = 0.52917721
    atc = np.random.random(size=(15,3))
    atq = np.random.random(size=15)
    gridc = np.random.random(size=(300,3))

    D = cdist(atc, gridc)
    Q = np.repeat(atq.reshape(-1,1), gridc.shape[0], axis=1)
    Vp = ( Q / (D * ang2au) ).sum(axis=0)

    Vf = getpot(atc, atq, gridc)

    print(np.allclose(Vp, Vf))
