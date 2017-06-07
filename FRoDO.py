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
    frlon = np.linspace((2*pi/f(2*frdim[1])), (2*pi)-(2*pi/(2*frdim[1])), num=frdim[1], dtype=np.double)
    frlat = np.linspace(-1+(1/(2*frdim[0])), 1-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)

    frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
    frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

    # Interpolate the magnetic field array to the output grid
    f0b = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), np.rot90(br[:,:,0]), kind='cubic')
    br0 = f0b(frlon, frlat)

    # Trace a uniform set of fieldlines for detection
    afl_r = np.zeros(frdim[0]*frdim[1])+1
    afl_ph, afl_th = np.meshgrid(frlon, frlat)
    afl_th = np.ndarray.flatten(afl_th)
    afl_ph = np.ndarray.flatten(afl_ph)

    afl_th = np.arcsin(-1*afl_th)+(pi/2)

    hlcy, fl = b.fieldlines_hlcy(afl_r, afl_th, afl_ph)
    
    # Begin by computing mappings of fieldline helicity and radial extent
    frhlcy[:,:] = 0
    frrext[:,:] = 0
    for ifl in np.arange(len(hlcy)):
        alfpr1 = fl[ifl][0][0]
        alfpr2 = fl[ifl][0][-1]
        alfpth1 = np.sin(fl[ifl][1][0] + pi/2)
        alfpth2 = np.sin(fl[ifl][1][-1] + pi/2)
        alfpph1 = fl[ifl][2][0]
        alfpph2 = fl[ifl][2][-1]
        
        alg1 = fl[ifl][0][0] < 1.2
        alg2 = fl[ifl][0][-1] < 1.2
        
        flrext = fl[ifl][0].max()
        
        why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
        whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
        why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
        whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
        if (len(whx1 == 1))&(len(why1) == 1)&(alg1):
            if abs(frhlcy[why1, whx1]) < abs(hlcy[ifl]):
                frhlcy[why1, whx1] = hlcy[ifl]
            if abs(frrext[why1, whx1]) < flrext:
                frrext[why1, whx1] = flrext
    if (len(whx2 == 1))&(len(why2) == 1)&(alg2):
        if abs(frhlcy[why2, whx2]) < abs(hlcy[ifl]):
            frhlcy[why2, whx2] = hlcy[ifl]
        if abs(frrext[why2, whx2]) < flrext:
            frrext[why2, whx2] = flrext

    # Begin the process of seeding flux rope footprints within this mapping
    havg = abs(frhlcy).mean()

    sthresh = (havg / ref_havg) * ref_sthresh
    ethresh = (havg / ref_havg) * ref_ethresh

    sthreshs[dcount] = sthresh
    ethreshs[dcount] = ethresh

    # Compute a set of core footprint regions
    sdfp = abs(frhlcy) > sthresh
    s = [[1,1,1],[1,1,1],[1,1,1]]
    sdfp_lab = (scipy.ndimage.label(sdfp, s))[0]

    for j in np.arange(sdfp_lab.max())+1:
        wharr = np.where(sdfp_lab == j)
        if len(wharr[0]) < 10:
            sdfp_lab[wharr] = 0
    sdfp_lab = (scipy.ndimage.label(sdfp_lab != 0, s))[0]

    # Compute a set of extent footprint regions
    expfp = abs(frhlcy) > ethresh
    exfp_lab = (scipy.ndimage.label(expfp, s))[0]

    for j in np.arange(expfp_lab.max())+1:
        wharr = np.where(expfp_lab == j)
            if len(wharr[0]) < 10:
                expfp_lab[wharr] = 0

    # Merge the core and extent regions
    fp_lab = np.zeros([len(frlat), len(frlon)], dtype='int32')
    for fri in np.arange(sdfp_lab.max())+1:
        fp_lab[np.where(expfp_lab == expfp_lab[np.where(sdfp_lab == fri)][0])] = 1

    # Separate adjacent regions with opposite sign of helicity
    frpos_lab = (scipy.ndimage.label((fp_lab!=0)*frhlcy>0, s))[0]
    for j in np.arange(frpos_lab.max())+1:
        wharr = np.where(frpos_lab == j)
        if len(wharr[0]) < 10:
            frpos_lab[wharr] = 0
    frpos_lab = (scipy.ndimage.label(frpos_lab != 0, s))[0]

    frneg_lab = (scipy.ndimage.label((fp_lab!=0)*frhlcy<0, s))[0]
    for j in np.arange(frneg_lab.max())+1:
        wharr = np.where(frneg_lab == j)
        if len(wharr[0]) < 10:
            frneg_lab[wharr] = 0
        frneg_lab = (scipy.ndimage.label(frneg_lab != 0, s))[0]

    frmap[:,:] = (frpos_lab+frneg_lab.max())*(frpos_lab!=0) + frneg_lab

    # Calculate footprint connectivity, and an index of flux rope fieldlines
    frcmap = np.zeros(frdim, dtype=np.int32)
    fr_fl = np.zeros(len(fl))
    for ifl in np.arange(len(hlcy)):
        alfpr1 = fl[ifl][0][0]
        alfpr2 = fl[ifl][0][-1]
        alfpth1 = np.sin(fl[ifl][1][0] + pi/2)
        alfpth2 = np.sin(fl[ifl][1][-1] + pi/2)
        alfpph1 = fl[ifl][2][0]
        alfpph2 = fl[ifl][2][-1]
        
        alg1 = fl[ifl][0][0] < 1.2
        alg2 = fl[ifl][0][-1] < 1.2
        
        why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
        whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
        why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
        whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
        
        ## Map out the pixel connectivity
        if (len(whx1 == 1))&(len(why1) == 1)&(frmap[why1,whx1] != 0):
            if (len(whx2 == 1))&(len(why2) == 1):
                frcmap[why1,whx1] = frmap[why2[0], whx2[0]]
                fr_fl[ifl] = frmap[why1, whx1]
        if (len(whx2 == 1))&(len(why2) == 1)&(frmap[why2,whx2] != 0):
            if (len(whx1 == 1))&(len(why1) == 1):
                frcmap[why2,whx2] = frmap[why1[0], whx1[0]]
                fr_fl[ifl] = frmap[why2, whx2]

    # Output appropriate arrays to disk
    outfile = netcdf.netcdf_file(outdir + 'fr_' + csfrm + '.nc', 'w')
    outfile.history = 'FRoDO flux rope data'
    outfile.createDimension('lat', frdim[0])
    outfile.createDimension('lon', frdim[1])
    
    out_frmap = outfile.createVariable('frmap', np.int16, ('lat', 'lon'))
    out_frmap[:] = frmap
    out_frmap.units = 'Flux rope footprint label (unitless)'
    
    out_frcmap = outfile.createVariable('frcmap', np.int16, ('lat', 'lon'))
    out_frcmap[:] = frcmap
    out_frcmap.units= 'Flux rope footprint connectivity label (unitless)'
    
    out_frhlcy = outfile.createVariable('frhlcy', np.double, ('lat', 'lon'))
    out_frhlcy[:] = frhlcy
    out_frhlcy.units = 'Field line cumulative helicity map (Mx cm^-2)'
    
    out_frrext = outfile.createVariable('frrext', np.double, ('lat', 'lon'))
    out_frrext[:]  = frrext
    out_frrext.units = 'FIeld line maximum radial extent (R_sun)'
    
    out_br0 = outfile.createVariable('br0', np.double, ('lat', 'lon'))
    out_br0[:] = br0
    out_br0.units = 'Radial magnetic flux density at 1.0 R_sun (G)'
    
    out_lat = outfile.createVariable('lat', np.float32, ('lat',))
    out_lat[:] = frlat
    out_lon = outfile.createVariable('lon', np.float32, ('lon',))
    out_lon[:] = frlon

    # Diagnostic readouts!
    time1 = datetime.now()
    if dcount == 0
        timedel = (time1 - time0)
    else:
        timedel = ((time1 - time0) + timedel) / 2
    timeeta = (nfrm - dcount) * timedel + time1
    print('Frame ' + '%05.f'%dcount + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta))

# Output any completed time-series variables
