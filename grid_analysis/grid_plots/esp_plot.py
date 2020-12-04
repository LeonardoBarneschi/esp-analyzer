#!/usr/bin/env python3
# coding: utf-8

"""
General Plot Parameters Dictionary
----------------------------------
- 'canvas_dimensions'   --> 'figsize' -- size of the plot (xy)
- 'x_title_dimension'   --> 'size' -- size of the x axis title
- 'y_title_dimension'   --> 'size' -- size of the x axis title
- 'z_title_dimension'   --> 'size' -- size of the x axis title
- 'output_quality'      --> 'dpi' -- density of pixels per inch
- 'labels_params'       --> 'rotation' -- label text angle
                        --> 'fontproperties' -- fonttype, dimensions, etc.
- 'ticks_params'        --> 'direction' -- tick facing inward or outward
                        --> 'width' -- tick width
                        --> 'length' -- tick length
                        --> 'color' -- tick color
- 'minor_ticks_params'  --> 'which' -- select 1 (x, y) or 2 axis (xy)
                        --> 'width' -- minor tick width
                        --> 'color' -- minor tick color
                        --> 'direction' -- minor tick direction

Contour Plot Parameters Dictionary
----------------------------------
- 'vmin' -- contour level minimum
- 'vmax' -- contour level maximum
- 'nlev' -- number of contour levels
- 'cmap' -- selected color map
- 'ccol' -- contour line color
- 'cwid' -- contour line width
- 'cfon' -- fontsize of the contour line level
- 'scol' -- scatter points color
- 'swid' -- scatter points size
"""

import os
import sys
import numpy as np
import pandas as pd
import csv
from scipy.interpolate import griddata
import matplotlib
import matplotlib.font_manager as fm
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter

general_params = {'canvas_dimensions'  : {'figsize':(6.0,4.0)},
                  'x_title_dimension'  : {'size':12},
                  'y_title_dimension'  : {'size':12},
                  'z_title_dimension'  : {'size':12},
                  'output_quality'     : {'dpi' : 600},
                  'labels_params'      : {'rotation' : 'horizontal','fontproperties' : fm.FontProperties(size=12)},
                  'ticks_params'       : {'direction' : 'in', 'axis' : 'both', 'width' : 1,'length' : 3,'color' : '#000000'},
                  'minor_ticks_params' : {'which' : 'minor','width' : 0.6,'color' : '#000000','direction' : 'in'}}

#contour_params = {'vmin': -5,
#                  'vmax':  5,
#contour_params = {'vmin': -4,
#                  'vmax':  4,
contour_params = {'vmin': -4,
                  'vmax':  0,
                  'nlev': 200,
                  'cmap':'RdBu',
                  #'cmap':'Spectral',
                  'cmap':'PiYG',
                  #'cmap':'RdYlBu',
                  'ccol': '#000000',
                  'cwid': 0.2,
                  'cfon': 8,
                  'scol': '#000000',
                  'swid': 30}

def myround(x, base=5):
    return base * round(float(x) / base)

def create_contour(x,y,z,resolution=100,contour_method='cubic'):

    resolution = str(resolution)+'j'
    X,Y = np.mgrid[min(x):max(x):complex(resolution), min(y):max(y):complex(resolution)]
    points = [[a,b] for a,b in zip(x,y)]
    Z = griddata(points, z, (X, Y), method=contour_method)
    return X,Y,Z

def add_vdw_contour(pts, ax, fig, png):

    ax.scatter(pts[:,0], pts[:,1], c='black',alpha=1, s=1)
    #ax.scatter(pts[:,0], pts[:,1], c='black',alpha=1, s=1)
    fig.savefig(png, bbox_inches='tight', **general_params['output_quality'])
    return

def plot_slice(df1, png, df2=False, tex=False):

    if tex:
        matplotlib.use('pgf')
        matplotlib.rcParams['pgf.texsystem'] = 'xelatex'
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.preamble'] = [ r'\usepackage{siunitx}',
                                                       r'\sisetup{detect-all}',
                                                       r'\usepackage{helvet}',
                                                       r'\usepackage{sansmath}',
                                                       r'\sansmath' ]
    else:
        matplotlib.use('TkAgg')

    A, B, C = create_contour(df1[1], df1[2], df1[3],resolution=100,contour_method='cubic')

    fig, ax = plt.subplots(**general_params['canvas_dimensions'])

    vmin, vmax, nlev, cmap, ccol, cwid, cfon, scol, swid = contour_params.values()

    # create dummy invisible image to use colormap desired on colorbar)
    img = plt.imshow(np.array([[vmin,vmax]]), cmap=plt.cm.get_cmap(cmap, nlev), aspect='auto')
    img.set_visible(False)

    cbar = plt.colorbar(orientation="vertical")
    cbar.ax.tick_params(direction='in')
    #tick_locator = ticker.MaxNLocator(nbins=nlev/2)
    #cbar.locator = tick_locator
    #cbar.update_ticks()

    ctr = ax.contourf(A, B, C, np.linspace(vmin, vmax, nlev), cmap=cmap, vmin=vmin, vmax=vmax, extend='both')


    try:
        xl, yl, zl = df2.min(axis=0).values
        xu, yu, zu = df2.max(axis=0).values
        for i in range(len(df2[0])):
            ax.scatter(df2[1][i], df2[2][i], color=scol, s=swid)
    except:
        pass

    if nlev <= 100:
        contours = ax.contour(A, B, C, np.linspace(vmin, vmax, nlev), colors=ccol,
                              linewidths=cwid)

        for line, lvl in zip(contours.collections, contours.levels):
            line.set_linestyle('-')

    # axis formtting 
    try:
        yu = myround(yu, 3)
        yl = myround(yl, 3)
        zu = myround(zu, 3)
        zl = myround(zl, 3)
        ax.set_xlim(yl, yu)
        ax.set_ylim(zl, zu)
    except:
        ax.set_xlim(-10, 8)
        ax.set_ylim(-3, 3)
    #ax.set_xlim(-14, 14)
    #ax.set_ylim(-14, 14)
    ax.minorticks_on()
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(**general_params['minor_ticks_params'])
    ax.tick_params(**general_params['ticks_params'])
    ax.set_yticklabels(ax.get_yticks(), **general_params['labels_params'])
    ax.set_xticklabels(ax.get_xticks(), **general_params['labels_params'])
    ax.yaxis.set_major_formatter(FormatStrFormatter('$%.0f$'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('$%.0f$'))
    ax.invert_yaxis()
    ax.invert_xaxis()

    # Title
    if tex==True:
        ax.set_xlabel('$\mathrm{X\hspace{0.01in}(\AA)}$')
        ax.set_ylabel('$\mathrm{Y\hspace{0.01in}(\AA)}$')
        cbar.set_label('$\mathrm{V_{Coul}\hspace{0.01in}(V)}$', labelpad=10)
    else:
        ax.set_xlabel('$\mathrm{X (\AA)}$')
        ax.set_ylabel('$\mathrm{Y (\AA)}$')
        cbar.set_label('$\mathrm{V_{Coul}\hspace{0.01in}(V)}$')
        plt.show()

    if nlev <= 20 and tex==False:
        plt.clabel(contours, inline=True, fontsize=cfon, fmt='%1.3f', use_clabeltext=True)

    fig.savefig(png, bbox_inches='tight', **general_params['output_quality'])

    return ax, fig

def parse_and_plot(esp, png, xyz=False, tex=False):

    df1 = pd.read_csv(esp, sep="\s+", skiprows=1, header=None)

    if xyz:
        df2 = pd.read_csv(xyz, sep="\s+", skiprows=1, header=None)
        ax, fig = plot_slice(df1, png, df2, tex)

    else:
        ax, fig = plot_slice(df1, png, tex=tex)

    return ax, fig, png

if __name__ == "__main__":
    pass
