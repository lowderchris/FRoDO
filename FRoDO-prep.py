# FRoDO-prep.py
# Subroutine to prepare input magnetic field data for flux rope detection
# In particular, the vector magnetic potential is computed in the DeVore gauge

# Import libraries
import os
import glob

import numpy as np
from scipy.io import netcdf

import compA

# Generate a list of files to search through
files = glob.glob('./dat/*.nc')
files.sort()

# Cycle through input data
for file in files:

    # Compute and store magnetic vector potential
    # Note that for now, for this naming scheme, just extract the timing string:
    # This will be fixed in future for more general data usage
    compA.compa(file[-11:-3])
