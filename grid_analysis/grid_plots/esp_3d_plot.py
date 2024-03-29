import os
import matplotlib
import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from scipy.interpolate import griddata
from scipy.spatial.distance import cdist
import itertools
import vtk

# Recorded script from Mayavi2
from numpy import array
mlab.options.offscreen = True

from mayavi.api import OffScreenEngine

def setview1():
    engine = OffScreenEngine()
    engine.start()
    #try:
    #    #engine = mayavi.engine
    #except NameError:
    #    from mayavi.api import Engine
    #    engine = Engine()
    #    engine.start()
    if len(engine.scenes) == 0:
        engine.new_scene()
    # -------------------------------------------
    scene = engine.scenes[0]
    scene.scene.background = (1.0, 1.0, 1.0)
    scene.scene.anti_aliasing_frames = 20
    scene.scene.parallel_projection = True
    scene.scene.camera.position = [40.575009454152, 5.59012378354684, -20.509715892269604]
    scene.scene.camera.focal_point = [0.36367527019261114, -1.7649274227504446, -0.46101479658482436]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.4233049466436812, -0.12981392512907605, -0.8966388721160671]
    scene.scene.camera.clipping_range = [30.556923337876515, 64.54772794093209]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    return scene

def setview2():
    if len(engine.scenes) == 0:
        engine.new_scene()
    # -------------------------------------------
    scene = engine.scenes[0]
    scene.scene.background = (1.0, 1.0, 1.0)
    scene.scene.camera.position = [-41.1699619564262, -14.827168967940311, 12.919429705884871]
    scene.scene.camera.focal_point = [0.36367527019261114, -1.7649274227504446, -0.46101479658482436]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.4233049466436812, -0.12981392512907605, -0.8966388721160671]
    scene.scene.camera.clipping_range = [30.556923337876515, 64.54772794093209]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()

def triarea(p1, p2, p3):

    a = np.linalg.norm(p1-p2)
    b = np.linalg.norm(p1-p3)
    c = np.linalg.norm(p2-p3)
    semi = (a + b + c) / 2
    area = semi*(semi-a)*(semi-b)*(semi-c) ** 0.5
    return area

def test(coords, p, png):

    D = cdist(coords, coords)
    M = np.zeros_like(D)
    for i in range(D.shape[0]):
        M[i] = D[i].argsort()

    triplets = []
    for i in range(M.shape[0]):
        triareas = []
        idx = M[i,0]
        combs = list(itertools.combinations(M[i,1:14], 2))
        triangles = [[ idx, comb[0], comb[1]] for comb in combs ]
        for t in triangles:
            p1, p2, p3 = coords[np.array(t).astype(int)]
            area = triarea(p1, p2, p3)
            triareas.append(area)

        triareas = np.array(triareas)
        triareas = np.argsort(triareas)
        triplets.extend(triangles[i] for i in triareas[-22:])
        #triplets.extend(triangles[i] for i in triareas[:15])
        #triplets.extend(triangles[i] for i in triareas[-22:-15])
        #triplets.extend(triangles[i] for i in triareas[-20:-13])


    fig = mlab.figure(size=(600, 600))
    scene = setview1()

    triplets = np.unique(triplets, axis=0).astype(int)
    print('Used {0} triplets.\n'.format(len(triplets)))
    m = mlab.triangular_mesh(coords[:,0], coords[:,1], coords[:,2], triplets, opacity=1, line_width=1, colormap="RdBu", scalars=p)
    mlab.savefig(png, magnification=4)

    #setview2()
    #print('Used {0} triplets.\n'.format(len(triplets)))
    #m = mlab.triangular_mesh(coords[:,0], coords[:,1], coords[:,2], triplets, opacity=1, line_width=1, colormap="RdBu", scalars=p)
    #mlab.savefig('/home/leonardo/surfaces/surface_back.png', magnification=10)

if __name__ == "__main__":

    path = os.getcwd()
    for fl in os.listdir(path):
        if fl.endswith('esp'):
            fl = inp
            data = np.loadtxt(inp, skiprows=1)
            coords = data[:,:3]
            p = data[:,-1]
            test(coords, p, fl[:-3]+'png')

    # inp = sys.argv[1]
    # coords = data[:,:3]
    # p = data[:,-1]
    # test(coords, p)
