# FRoDO-erupt.py
# Flux Rope Detection and Organization
# Subroutine for the detection, tracking, and linkage of erupting flux ropes with the full database of ropes generated from FRoDO.py

# Import libraries
import b_sim_netcdf
import os
import glob
import datetime
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
bdatprefix = config['paths']['bdatprefix']
adatprefix = config['paths']['adatprefix']
outdir = config['paths']['outdir']

frdim = np.array([np.int(config['array']['nlat']), np.int(config['array']['nlon'])])
maxlat = np.double(config['array']['maxlat'])

ref_sthresh = np.double(config['thresholds']['ref_sthresh'])
ref_ethresh = np.double(config['thresholds']['ref_ethresh'])
ref_bavg = np.double(config['thresholds']['ref_bavg'])
ref_havg = np.double(config['thresholds']['ref_havg'])

# Generate a list of files to search through
bfrm_list = glob.glob(datdir + bdatprefix + '*.nc')
bfrm_list.sort()
afrm_list = glob.glob(datdir + adatprefix + '*.nc')
afrm_list.sort()
nfrm = len(bfrm_list)

# Initialize arrays for storage
tarr = []

fr_elab = np.array([])
fr_efpt = np.array([])
fr_etarr = np.array([])

phfl0 = np.zeros(frdim[0]*frdim[1],dtype=np.double)
thfl0 = np.zeros(frdim[0]*frdim[1],dtype=np.double)
rfl0 = np.zeros(frdim[0]*frdim[1],dtype=np.double)
phfl1 = np.zeros(frdim[0]*frdim[1],dtype=np.double)
thfl1 = np.zeros(frdim[0]*frdim[1],dtype=np.double)
rfl1 = np.zeros(frdim[0]*frdim[1],dtype=np.double)

# Begin cycling through these frames
dcount = 0
for cfrm in bfrm_list:

    # Define some timing
    time0 = datetime.datetime.now()
    csfrm = cfrm[-11:-3]

    # Read original data into memory
    b = b_sim_netcdf.SphB_sim(datdir + bdatprefix + csfrm + '.nc', datdir + adatprefix + csfrm + '.nc', 128,128,128)
    d = netcdf.netcdf_file(datdir + bdatprefix + csfrm + '.nc', 'r')

    r = d.variables['r'][:].copy()
    th = d.variables['th'][:].copy()
    ph = d.variables['ph'][:].copy()

    br = d.variables['br'][:,:,:].copy()
    bth = d.variables['bth'][:,:,:].copy()
    bph = d.variables['bph'][:,:,:].copy()
    cdate = d.date
    ctarr = datetime.datetime.strptime(bytes.decode(cdate),"%Y%b%d_%H%M")
    cjuld = np.double(sunpy.time.julian_day(ctarr))
    tarr.append(ctarr)

    # Close netCDF files
    d.close()

    # Read processed data into memory
    d = netcdf.netcdf_file(outdir + 'fr-' + csfrm + '.nc', 'r')
    frhlcy = d.variables['frhlcy'][:,:].copy()
    frmap = d.variables['frmap'][:,:].copy()
    d.close()

    # Read the previous timeframe into memory for comparison
    if dcount != 0:
        d = netcdf.netcdf_file(outdir + 'fr-' + csfrm1 + '.nc', 'r')
        frmap0 = d.variables['frmap'][:,:].copy()
        d.close()

    # Define some coordinate information
    lons, lats = np.meshgrid(ph*360./(2*pi), th*360./(2*pi)-90.)
    frlon = np.linspace((2*pi/(2*frdim[1])), (2*pi)-(2*pi/(2*frdim[1])), num=frdim[1], dtype=np.double)
    frlat = np.linspace(-np.sin(maxlat*np.pi/180.)+(1/(2*frdim[0])), np.sin(maxlat*np.pi/180.)-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)
    frlat_pol = np.linspace(-1+(1/(2*frdim[0])), 1-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)

    frlat_wcell = abs(frlat[0] - frlat[1])
    frlon_wcell = abs(frlon[0] - frlon[1])

    frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
    frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

    # Interpolate some of this data to the output grid
    bh_intarr = np.rot90((bth[:,:,-1]**2 + bph[:,:,-1]**2)**(0.5))
    f0bh = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), bh_intarr, kind='cubic')
    bht = f0bh(frlon, frlat)
    bht_pol = f0bh(frlon, frlat_pol)

    # Trace a set of fieldlines down from the simulation upper boundary
    afl_r = np.zeros(frdim[0]*frdim[1])+2.5
    afl_ph, afl_th = np.meshgrid(frlon, frlat)
    afl_th = np.ndarray.flatten(afl_th)
    afl_ph = np.ndarray.flatten(afl_ph)

    afl_th = np.arcsin(-1*afl_th)+(pi/2)

    hlcy, fl = b.fieldlines_hlcy(afl_r, afl_th, afl_ph)

    ## Generate an array of field-line helicity
    hlcy_tb = np.zeros([frdim[0], frdim[1]], dtype='float64')
    for fl_ind in np.arange(len(hlcy)):
        alfpr1 = fl[fl_ind][0][0]
        alfpr2 = fl[fl_ind][0][-1]
        alfpth1 = np.sin(fl[fl_ind][1][0] + pi/2)
        alfpth2 = np.sin(fl[fl_ind][1][-1] + pi/2)
        alfpph1 = fl[fl_ind][2][0]
        alfpph2 = fl[fl_ind][2][-1]

        alg1 = fl[fl_ind][0][0] > 2.4
        alg2 = fl[fl_ind][0][-1] > 2.4

        why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
        whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
        why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
        whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
        if (len(whx1 == 1))&(len(why1) == 1)&(alg1):
            if abs(hlcy_tb[why1, whx1]) < abs(hlcy[fl_ind]):
                hlcy_tb[why1, whx1] = hlcy[fl_ind]
        if (len(whx2 == 1))&(len(why2) == 1)&(alg2):
            if abs(hlcy_tb[why2, whx2]) < abs(hlcy[fl_ind]):
                hlcy_tb[why2, whx2] = hlcy[fl_ind]

    # Compute a few parameters for thresholding
    bavg = abs(bht).mean()
    havg = abs(frhlcy).mean()
    shthresh = (havg / 0.26648107322354925) * 1.00
    ehthresh = (havg / 0.26648107322354925) * 0.70
    sbthresh = (bavg / 0.027634660528377285) * 0.30
    ebthresh = (bavg / 0.027634660528377285) * 0.10

    # Sort out regions with a reasonable overlap of horizontal magnetic field and helicity
    regb = scipy.ndimage.label(bht > ebthresh)[0]

    # Filter out small regions
    for r in np.arange(regb.max())+1:
        wregb = np.where(regb ==  r)
        if len(wregb[0]) < 100:
            regb[wregb] = 0
        if abs(hlcy_tb[wregb]).max() < ehthresh:
            regb[wregb] = 0

    # Relabel regions
    regb = scipy.ndimage.label(regb!=0)[0]

    # Grab a mapping of fieldline start/end points
    for fli in np.arange(len(fl)):
        phfl0[fli] = fl[fli][2][0]
        phfl1[fli] = fl[fli][2][-1]
        thfl0[fli] = fl[fli][1][0]
        thfl1[fli] = fl[fli][1][-1]
        rfl0[fli] = fl[fli][0][0]
        rfl1[fli] = fl[fli][0][-1]

    # Define a list of fieldlines that fall within this region
    regls0 = np.array([], dtype=np.int)
    regls1 = np.array([], dtype=np.int)
    regfl0 = np.array([], dtype=np.int)
    regfl1 = np.array([], dtype=np.int)

    # Alternatively, define a set of fieldlines just by grabbing the line intersecting closest to the cell center
    # Find the corresponding fieldlines contained within, and trace back to origins 
    for r in np.arange(regb.max())+1:
        wregb = np.where(regb == r)
        for rci in np.arange(len(wregb[0])):
            rlat = wregb[0][rci]
            rlon = wregb[1][rci]
            fltb0 = np.where(rfl0 > 2.4)
            fltb1 = np.where(rfl1 > 2.4)

            londff0 = abs(phfl0[fltb0] - frlon[rlon])
            latdff0 = abs(thfl0[fltb0] - (np.arcsin(-1*frlat[rlat])+(pi/2)))
            londff1 = abs(phfl1[fltb1] - frlon[rlon])
            latdff1 = abs(thfl1[fltb1] - (np.arcsin(-1*frlat[rlat])+(pi/2)))
            # Note that this 'distance' approximation computed here is only for comparison. A true distance calculation could be substituted here, when consolidating codes together.
            fldist0 = np.sqrt(londff0**2 + latdff0**2)
            fldist1 = np.sqrt(londff1**2 + latdff1**2)
            wmd0 = np.where(fldist0 == np.min(fldist0))[0]
            wmd1 = np.where(fldist1 == np.min(fldist1))[0]
            if ((londff0[wmd0] < frlon_wcell) & (latdff0[wmd0] < frlat_wcell)):
                regls0 = np.append(regls0,r)
                regfl0 = np.append(regfl0, np.where(phfl0 == phfl0[fltb0][wmd0])[0][0])
                #regfl0 = np.append(regfl0, np.where(phfl0 == phfl0[fltb0][where(fldist0 == mindist0)])[0][0])
                # Go ahead here and assign this value to the cell below? Will need to compute coordinates?
            if ((londff1[wmd1] < frlon_wcell) & (latdff1[wmd1] < frlat_wcell)):
                regls1 = np.append(regls1,r)
                regfl1 = np.append(regfl1, np.where(phfl1 == phfl1[fltb1][wmd1])[0][0])
                #regfl1 = np.append(regfl1, np.where(phfl1 == phfl1[fltb1][where(fldist1 == mindist1)])[0][0])

    # To find corresponding footpoint locations, simply invert the indicies and check for radius value
    # Create a linked map at 1R_sun to contain footprints of erupting flux ropes

    regbb = np.zeros([frdim[0], frdim[1]], dtype=np.double)
    regc = 0
    for fl_ind in regfl0:
        alfpr1 = fl[fl_ind][0][0]
        alfpr2 = fl[fl_ind][0][-1]
        alfpth1 = np.sin(fl[fl_ind][1][0] + pi/2)
        alfpth2 = np.sin(fl[fl_ind][1][-1] + pi/2)
        alfpph1 = fl[fl_ind][2][0]
        alfpph2 = fl[fl_ind][2][-1]

        alg1 = fl[fl_ind][0][0] < 1.2
        alg2 = fl[fl_ind][0][-1] < 1.2

        if alg1:
            why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
            whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
            if (len(whx1) == 1)&(len(why1) == 1)&(alg1):
                regbb[why1, whx1] = regls0[regc]
        if alg2:
            why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
            whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
            if (len(whx2) == 1)&(len(why2) == 1)&(alg2):
                regbb[why2, whx2] = regls0[regc]
        regc = regc + 1
    regc = 0
    for fl_ind in regfl1:
        alfpr1 = fl[fl_ind][0][0]
        alfpr2 = fl[fl_ind][0][-1]
        alfpth1 = np.sin(fl[fl_ind][1][0] + pi/2)
        alfpth2 = np.sin(fl[fl_ind][1][-1] + pi/2)
        alfpph1 = fl[fl_ind][2][0]
        alfpph2 = fl[fl_ind][2][-1]

        alg1 = fl[fl_ind][0][0] < 1.2
        alg2 = fl[fl_ind][0][-1] < 1.2

        if alg1:
            why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
            whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
            if (len(whx1) == 1)&(len(why1) == 1)&(alg1):
                regbb[why1, whx1] = regls1[regc]
        if alg2:
            why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
            whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
            if (len(whx2) == 1)&(len(why2) == 1)&(alg2):
                regbb[why2, whx2] = regls1[regc]
        regc = regc + 1

    # Compute a pixel shift amount to account for differential rotation
    if dcount == 0:
        dtime = 1
    else:
        dtime0 = tarr[dcount] - tarr[dcount - 1]
        dtime = dtime0.days + (dtime0.seconds / 86400.)
    drot = diff_rot(dtime * u.day, np.arcsin(frlat)*180/pi * u.deg, rot_type='snodgrass')
    drot = drot - drot.max()
    drot = np.array(np.around((drot/u.deg * 180/pi) * (frlon[1]-frlon[0]))).astype(np.int)

    # Begin comparison with detected flux rope footprints 
    # Compute arrays storing the locations of flux rope footprints
    fpreg = np.array([])
    fpwhr = list([])
    if dcount == 0:
        for fr in np.unique(frmap):
            if fr != 0:
                fpreg = np.append(fpreg,fr)
                fpwhr.append(np.where(frmap == fr))
    else:
        for r in np.arange(frdim[0]):
            frmap0[r,:] = np.roll(frmap0[r,:], drot[r])
        for fr in np.unique(frmap0):
            if fr != 0:
                fpreg = np.append(fpreg,fr)
                fpwhr.append(np.where(frmap == fr))

    # Run through this set and compare for overlap
    for fr in np.arange(regb.max())+1:
        frwhr = np.where(regb == fr)
        frwhrbb = np.where(regbb == fr)
        fpmaxarea = 0.
        fpmaxarea_reg = 0.
        # Run through and compare overlap with full flux rope dataset
        for ifpfr in np.arange(len(fpreg)):
            xlap = np.in1d(fpwhr[ifpfr][1], frwhrbb[1])
            ylap = np.in1d(fpwhr[ifpfr][0], frwhrbb[0])
            plap = np.where(xlap & ylap)[0]
            if (len(plap) != 0) & (len(np.where(np.in1d(fr_elab, int(fpreg[ifpfr])))[0]) == 0):
                fr_elab = np.append(fr_elab, int(fpreg[ifpfr]))
                fr_etarr = np.append(fr_etarr, csfrm)
                if len(fpwhr[ifpfr][1]) > fpmaxarea:
                    fpmaxarea = len(fpwhr[ifpfr][1])
                    fpmaxarea_reg = int(fpreg[ifpfr])
        if (fpmaxarea_reg != 0):
            fr_efpt = np.append(fr_efpt, fpmaxarea_reg)

    # Diagnostic readouts!
    time1 = datetime.datetime.now()
    if dcount == 0:
        timedel = (time1 - time0)
    else:
        timedel = ((time1 - time0) + timedel) / 2
    timeeta = (nfrm - (dcount+1)) * timedel + time1
    print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta))

    # Save the current timestamp
    csfrm1 = csfrm

    # Advance the timing index
    dcount = dcount + 1

# Save this data to file
outfile = open(outdir + 'hist/fr-elab.pkl', 'wb')
pickle.dump(fr_elab, outfile)
outfile.close()

outfile = open(outdir + 'hist/fr-efpt.pkl', 'wb')
pickle.dump(fr_efpt, outfile)
outfile.close()

outfile = open(outdir + 'hist/fr-etarr.pkl', 'wb')
pickle.dump(fr_etarr, outfile)
outfile.close()
