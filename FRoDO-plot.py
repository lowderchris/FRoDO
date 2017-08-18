# FRoDO-plot.py
# Set of routines to plot more advanced plotting outputs
'''
This set of modules works alongside generated data from the FRoDO code.
Specify plotting parameters within the config.cfg file.
For all additional details, consult the aptly named README file.
'''

## To run headlessly in the background,
## xvfb-run --server-args="-screen 0 1024x768x24" python FRoDO-plot.py

# Import libraries
import os
import glob
import re
from datetime import datetime
import matplotlib
from matplotlib import cm, colors
import pickle
from scipy.io import netcdf
import numpy as np
from numpy import pi
import palettable
import b_sim_netcdf

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

datdir = config['paths']['datdir']
bdatprefix = config['paths']['bdatprefix']
adatprefix = config['paths']['adatprefix']
outdir = config['paths']['outdir']

frmdir = config['paths']['frmdir']
movname = config['paths']['movname']

flalph = config['plot3d']['flalph']
pltanno = config['plot3d']['pltanno']
pltlght = config['plot3d']['pltlght']
vph0 = config['plot3d']['vph0']
vph1 = config['plot3d']['vph1']
vth0 = config['plot3d']['vth0']
vth1 = config['plot3d']['vth1']
vr0 = config['plot3d']['vr0']
vr1 = config['plot3d']['vr1']
frmrate = config['plot3d']['frmrate']

def plot3d():
    '''
    Uses specified plotting parameters to animate a three-dimension visualization of magnetic flux rope evolution.
    '''

    # Import 3D plotting libraries
    from mayavi import mlab

    # Generate a list of files for animating
    files = glob.glob(outdir + 'fr-'+'*.nc')
    files.sort()

    # Read radial extent and duration filtered index
    infile = open(outdir + '/hist/fr-frg.pkl', 'rb')
    fr_frg = pickle.load(infile)
    infile.close()
    fr_frg = fr_frg.astype(np.int)

    # Iterate through these files
    for file in files:
        csfrm = re.split('fr-|\.', file)[1]

        # Read required data
        b = b_sim_netcdf.SphB_sim(datdir+bdatprefix+csfrm+'.nc', datdir+adatprefix+csfrm+'.nc', 128, 128, 128)
        f = netcdf.netcdf_file(outdir + 'fr-'+csfrm+'.nc', 'r')
        frmap = f.variables['frmap'][:,:].copy()
        br0 = f.variables['br0'][:,:].copy()
        frhlcy = f.variables['frhlcy'][:,:].copy()
        lat = f.variables['lat'][:].copy()
        lon = f.variables['lon'][:].copy()
        f.close()

        # Create array of filtered flux rope footprints to trace

        # Re-trace these fieldlines

        # Setup the scene and plot fieldlines

        # Animation additional frames if required

def plot2d():
    '''
    Flattens data to a two-dimensional representation for a latitude-longitude projection plot of flux rope fieldlines.
    '''

