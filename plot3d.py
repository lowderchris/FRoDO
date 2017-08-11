# plot3d.py
# Set of routines to plot FRoDO output data in
# the third dimension.

## To run headlessly in the background,
## xvfb-run --server-args="-screen 0 1024x768x24" python plot3d.py

# Import libraries
import os
import glob
from datetime import datetime
import matplotlib
from matplotlib import cm, colors
import pickle
from scipy.io import netcdf
import numpy as np
from numpy import pi
import palettable
import b_sim_netcdf
from mayavi import mlab

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

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

    
