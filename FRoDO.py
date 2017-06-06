# FRoDO.py
# Flux Rope Detection and Observation

# Import libraries
import b_sim_netcdf
import os
import glob
from datetime import datetime
import cPickle as pickle

import numpy as np
from numpy import pi
import scipy
from scipy.io import b_sim_netcdf

import sunpy.time

# Read a few user parameters from file
# Include here data/output directories
# For now, code these in place here
frdim = [180,360]               # Output array dimensions
sfrm = 904                      # Starting frame
efrm = 904                      # Ending frame
dfrm = 1                        # Frame step
ref_sthresh = 1.00              # Reference seed threshold
ref_ethresh = 0.70              # Reference extent threshold
ref_bavg = 4.5508878954543137	# Reference mean B
ref_havg = 0.26648107322354925	# Reference mean helicity

datdir = './dat/'
labstr = 's{0:.2f}e{1:.2f}'.format(ref_sthresh,ref_ethresh)
outdir = datdir + 'fr' + labstr + '/'

os.system("mkdir" + outdir)

# Work through all of the nasty date calculations
frm_list = np.arange(sfrm,efrm+1,step=dfrm)
nfrm = len(frm_list)

# Define any other variables needed
dcount = 0
frmap = np.zeros(frdim, dtype=np.int16)
frcmap = np.zeros(frdim, dtype=np.in16)
frhlcy = np.zeros(frdim, dtype=np.double)
frrext = np.zeros(frdim, dtype=np.double)

br0 = np.zeros(frdim, dtype=np.double)

tarr = []

sthreshs = np.zeros(nfrm, dtype=np.double)
ethreshs = np.zeros(nfrm, dtype=np.double)

# Begin cycling through time frames
for cfrm in frm_list:

    # Define some timing
    time0 = datetime.now()
    csfrm = '%05.f'%cfrm

    # Read data into memory
    b = b_sim_netcdf.SphB_sim(datdir + 'b_' + csfrm + '.nc', 128,128,128)
    d = netcdf.netcdf_file(datdir + 'b_' + csfrm + '.nc', 'r')

    r = d.variables['r'][:].copy()
    th = d.variables['th'][:].copy()
    ph = d.variables['ph'][:].copy()

    ar = d.variables['ar'][:,:,:].copy()
    ath = d.variables['ath'][:,:,:].copy()
    aph = d.variables['aph'][:,:,:].copy()
    br = d.variables['br'][:,:,:].copy()
    bth = d.variables['bth'][:,:,:].copy()
    bph = d.variables['bph'][:,:,:].copy()
    jr = d.variables['jr'][:,:,:].copy()
    jth = d.variables['jth'][:,:,:].copy()
    jph = d.variables['jph'][:,:,:].copy()
    cdate = d.date
    ctarr = datetime.strptime(cdate,"%Y%b%d_%H%M")
    cjuld = np.double(sunpy.time.julian_day(ctarr))
    tarr.append(ctarr)

    # Close netCDF files
    d.close()

    # Define some coordinate information
    lons, lats = np.meshgrid(ph*360./(2*pi), th*360./(2*pi)-90.)
    frlon = np.linspace((2*pi/720.), (2*pi)-(2*pi/720.), num=360, dtype=np.double)
    frlat = np.linspace(-1+(1/180.), 1-(1/180.), num=180, dtype=np.double)

    frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
    frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

    # Interpolate the magnetic field array to the output grid
    f0b = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), np.rot90(br[:,:,0]), kind='cubic')
    br0 = f0b(frlon, frlat)

    # Trace a uniform set of fieldlines for detection
    afl_r = np.zeros(180*360)+1
    afl_ph, afl_th = np.meshgrid(frlon, frlat)
    afl_th = np.ndarray.flatten(afl_th)
    afl_ph = np.ndarray.flatten(afl_ph)

    afl_th = np.arcsin(-1*afl_th)+(pi/2)

    hlcy, jpar, fl = b.fieldlines_hlcy_jpar(afl_r, afl_th, afl_ph)

	# Here there be flux rope detection

	# Output appropriate arrays to disk

    # Diagnostic readouts!
    time1 = datetime.now()
    if dcount == 0
        timedel = (time1 - time0)
    else:
        timedel = ((time1 - time0) + timedel) / 2
    timeeta = (nfrm - dcount) * timedel + time1
    print('Frame ' + '%05.f'%dcount + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta))

# Output any completed time-series variables
