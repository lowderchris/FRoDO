# FRoDO-prep.py
# Subroutine to prepare input magnetic field data for flux rope detection
# In particular, the vector magnetic potential is computed in the DeVore gauge

# Import libraries
import numpy as np
from scipy.io import netcdf

import compA

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

sfrm = np.int(config['times']['sfrm'])
efrm = np.int(config['times']['efrm'])
dfrm = np.int(config['times']['dfrm'])

# Generate a list of files to search through
frm_list = np.arange(sfrm,efrm+1,step=dfrm)
nfrm = len(frm_list)

# Cycle through input data
dcount = 0
for cfrm in frm_list:

    # Compute and store magnetic vector potential
    compA.compa(cfrm)
