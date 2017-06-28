# FRoDO-erupt.py
# Flux Rope Detection and Organization
# Subroutine for the detection, tracking, and linkage of erupting flux ropes with the full database of ropes generated from FRoDO.py

# Import libraries
import b_sim_netcdf
import os
import glob
from datetime import datetime
#import cPickle as pickle
import pickle

import numpy as np
from numpy import pi
import scipy
import scipy.stats
from scipy.io import netcdf

import sunpy.time
from sunpy.physics.differential_rotation import diff_rot
import astropy.units as u

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

datdir = config['paths']['datdir']
outdir = config['paths']['outdir']

sfrm = np.int(config['times']['sfrm'])
efrm = np.int(config['times']['efrm'])
dfrm = np.int(config['times']['dfrm'])

frdim = np.array([np.int(config['array']['nlat']), np.int(config['array']['nlon'])])

ref_sthresh = np.double(config['thresholds']['ref_sthresh'])
ref_ethresh = np.double(config['thresholds']['ref_ethresh'])
ref_bavg = np.double(config['thresholds']['ref_bavg'])
ref_havg = np.double(config['thresholds']['ref_havg'])

# Generate a list of files to search through
frm_list = np.arange(sfrm,efrm+1,step=dfrm)
nfrm = len(frm_list)

# Begin cycling through these frames
dcount = 0
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

    br = d.variables['br'][:,:,:].copy()
    bth = d.variables['bth'][:,:,:].copy()
    bph = d.variables['bph'][:,:,:].copy()
    cdate = d.date
    ctarr = datetime.strptime(bytes.decode(cdate),"%Y%b%d_%H%M")
    cjuld = np.double(sunpy.time.julian_day(ctarr))
    tarr.append(ctarr)

    # Close netCDF files
    d.close()

    # Define some coordinate information
    lons, lats = np.meshgrid(ph*360./(2*pi), th*360./(2*pi)-90.)
    frlon = np.linspace((2*pi/(2*frdim[1])), (2*pi)-(2*pi/(2*frdim[1])), num=frdim[1], dtype=np.double)
    frlat = np.linspace(-1+(1/(2*frdim[0])), 1-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)

    frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
    frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

    # Interpolate some of this data to the output grid
    bh_intarr = np.rot90((bth[:,:,-1]**2 + bph[:,:,-1]**2)**(0.5))
    f0bh = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), bh_intarr, kind='cubic')
    bht = f0bh(fplon, fplat)

    # Trace a set of fieldlines down from the simulation upper boundary
    afl_r = np.zeros(frdim[0]*frdim[1])+2.5
    afl_ph, afl_th = np.meshgrid(frlon, frlat)
    afl_th = np.ndarray.flatten(afl_th)
    afl_ph = np.ndarray.flatten(afl_ph)

    afl_th = np.arcsin(-1*afl_th)+(pi/2)

    hlcy, fl = b.fieldlines_hlcy(afl_r, afl_th, afl_ph)
