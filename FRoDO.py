# FRoDO.py
# Flux Rope Detection and Organization

# Import libraries
import b_sim_netcdf
import os
import glob
from datetime import datetime
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

sfrm = np.int(config['times']['sfrm'])
efrm = np.int(config['times']['efrm'])
dfrm = np.int(config['times']['dfrm'])

frdim = np.array([np.int(config['array']['nlat']), np.int(config['array']['nlon'])])

ref_sthresh = np.double(config['thresholds']['ref_sthresh'])
ref_ethresh = np.double(config['thresholds']['ref_ethresh'])
ref_bavg = np.double(config['thresholds']['ref_bavg'])
ref_havg = np.double(config['thresholds']['ref_havg'])

# Create output directories if needed
os.system("mkdir " + outdir)
os.system("mkdir " + outdir + 'hist/')

# Remove any existing files
os.system("rm " + outdir + "*")
os.system("rm " + outdir + "hist/*")

# Work through all of the nasty date calculations
frm_list = np.arange(sfrm,efrm+1,step=dfrm)
nfrm = len(frm_list)

# Define arrays to store surface maps
frmap = np.zeros(frdim, dtype=np.int16)
frcmap = np.zeros(frdim, dtype=np.int16)
frhlcy = np.zeros(frdim, dtype=np.double)
frrext = np.zeros(frdim, dtype=np.double)
br0 = np.zeros(frdim, dtype=np.double)

# Define an empty array of times
tarr = []

# Define a pixel area array (note the use of a sine latitude grid here)
rsun = 6.957e10
pix_area = 4. * np.pi * rsun**2 / (frdim[0] * frdim[1])

# Define arrays for the storage of flux rope time histories
fr_area = np.array([np.nan])
fr_time = np.array([np.nan])
fr_dur = np.array([np.nan])
fr_mhlcy = np.array([np.nan])
fr_nhlcy = np.array([np.nan])
fr_sflux = np.array([np.nan])
fr_uflux = np.array([np.nan])
fr_mlat = np.array([np.nan])
fr_rext = np.array([np.nan])
fr_mrext = np.array([np.nan])

frh_area = list([np.nan])
frh_time = list([np.nan])
frh_mhlcy = list([np.nan])
frh_nhlcy = list([np.nan])
frh_sflux = list([np.nan])
frh_uflux = list([np.nan])
frh_rext = list([np.nan])
frh_mrext = list([np.nan])

# Define arrays to store threshold levels
sthreshs = np.zeros(nfrm, dtype=np.double)
ethreshs = np.zeros(nfrm, dtype=np.double)

# Define some loop counting information
dcount = 0

# Begin cycling through time frames
for cfrm in frm_list:

    # Define some timing
    time0 = datetime.now()
    csfrm = '%05.f'%cfrm

    # Read magnetic field data into memory
    b = b_sim_netcdf.SphB_sim(datdir + bdatprefix + csfrm + '.nc', datdir + adatprefix + csfrm + '.nc', 128,128,128)
    d = netcdf.netcdf_file(datdir + bdatprefix + csfrm + '.nc', 'r')

    r = d.variables['r'][:].copy()
    th = d.variables['th'][:].copy()
    ph = d.variables['ph'][:].copy()

    ar = d.variables['ar'][:,:,:].copy()
    ath = d.variables['ath'][:,:,:].copy()
    aph = d.variables['aph'][:,:,:].copy()
    br = d.variables['br'][:,:,:].copy()
    bth = d.variables['bth'][:,:,:].copy()
    bph = d.variables['bph'][:,:,:].copy()
    cdate = d.date
    ctarr = datetime.strptime(bytes.decode(cdate),"%Y%b%d_%H%M")
    cjuld = np.double(sunpy.time.julian_day(ctarr))
    tarr.append(ctarr)

    d.close()

    # Read vector magnetic potential data into memory
    d = netcdf.netcdf_file(datdir + adatprefix + csfrm + '.nc', 'r')

    ar = d.variables['ar'][:,:,:].copy()
    ath = d.variables['ath'][:,:,:].copy()
    aph = d.variables['aph'][:,:,:].copy()

    d.close()

    # Define some coordinate information
    lons, lats = np.meshgrid(ph*360./(2*pi), th*360./(2*pi)-90.)
    frlon = np.linspace((2*pi/(2*frdim[1])), (2*pi)-(2*pi/(2*frdim[1])), num=frdim[1], dtype=np.double)
    frlat = np.linspace(-1+(1/(2*frdim[0])), 1-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)

    frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
    frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

    # Interpolate the magnetic field array to the output grid
    f0b = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), np.rot90(br[:,:,0]), kind='cubic')
    br0 = f0b(frlon, frlat)

    # Trace a uniform set of fieldlines for detection
    afl_r = np.zeros(frdim[0]*frdim[1])+1.0
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
    expfp_lab = (scipy.ndimage.label(expfp, s))[0]

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

    # Compute a pixel shift amount to account for differential rotation
    if dcount == 0:
        drot = diff_rot(1 * u.day, np.arcsin(frlat)*180/pi * u.deg, rot_type='snodgrass')
        drot = drot - drot.max()
        drot = np.array(np.around((drot/u.deg * 180/pi) * (frlon[1]-frlon[0]))).astype(np.int)

    frnlab = np.zeros(frdim, dtype=np.int16)
    # Compare flux rope maps with prior times to label histories
    if dcount != 0:
        frlab = np.copy(frmap)
        for r in np.arange(frdim[0]):
            frmap0[r,:] = np.roll(frmap0[r,:], drot[r])
            frhlcy0[r,:] = np.roll(frhlcy0[r,:], drot[r])
        # Label and iterate through flux rope footprint regions
        #for fr in np.arange(frlab.max())+1:
        for fr in np.unique(frlab):
            if fr != 0:
                frwhr = np.where(frlab == fr)
                frarr = frmap0[frwhr]
                hlcyarr0 = frhlcy0[frwhr]
                hlcyarr1 = frhlcy[frwhr]
                hlcycheck = (np.sign(hlcyarr0.mean()) * np.sign(hlcyarr1.mean())) == 1
                frsts = scipy.stats.mode(frarr)
                if (float(frsts[1][0]) / len(frarr) > 0.5) & (hlcycheck) & (frsts[0][0] != 0):
                    frnlab[frwhr] = frsts[0][0]
                else:
                    frnlab[frwhr] = regmx + 1
                    regmx = regmx + 1
        if frnlab.max() != 0 : regmx = frnlab.max()
        frmap = np.copy(frnlab)
    else:
        regmx = frmap.max()

    # Record some statistics
    for frcur in np.unique(frmap[frmap!=0]):
        # Output location information
        frwhr = np.where(frmap == frcur)
        if len(fr_area) <= frcur:
            fr_area = np.append(fr_area, len(frwhr[0]) * pix_area)
            fr_time = np.append(fr_time, dcount)
            fr_mlat = np.append(fr_mlat, np.mean(frlat[frwhr[0]]))
            fr_dur = np.append(fr_dur, 1)
            fr_mhlcy = np.append(fr_mhlcy, (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean())
            fr_nhlcy = np.append(fr_nhlcy, (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum())
            fr_sflux = np.append(fr_sflux, (br0[frwhr] * pix_area).sum())
            fr_uflux = np.append(fr_uflux, abs(br0[frwhr] * pix_area).sum())
            fr_rext = np.append(fr_rext, frrext[frwhr].mean())
            fr_mrext = np.append(fr_mrext, frrext[frwhr].max())

            frh_area.append(np.array([len(frwhr[0]) * pix_area]))
            frh_time.append(np.array([dcount]))
            frh_mhlcy.append(np.array([(hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean()]))
            frh_nhlcy.append(np.array([(hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum()]))
            frh_sflux.append(np.array([(br0[frwhr] * pix_area).sum()]))
            frh_uflux.append(np.array([abs(br0[frwhr] * pix_area).sum()]))
            frh_rext.append(np.array([frrext[frwhr].mean()]))
            frh_mrext.append(np.array([frrext[frwhr].max()]))
        else:
            cflux = abs(br0[frwhr] * pix_area).sum()
            fr_dur[frcur] = fr_dur[frcur] + 1

            frh_area[frcur] = np.append(frh_area[frcur], len(frwhr[0]) * pix_area)
            frh_time[frcur] = np.append(frh_time[frcur], dcount)
            frh_mhlcy[frcur] = np.append(frh_mhlcy[frcur], (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean())
            frh_nhlcy[frcur] = np.append(frh_nhlcy[frcur], (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum())
            frh_sflux[frcur] = np.append(frh_sflux[frcur], (br0[frwhr] * pix_area).sum())
            frh_uflux[frcur] = np.append(frh_uflux[frcur], abs(br0[frwhr] * pix_area).sum())
            frh_rext[frcur] = np.append(frh_rext[frcur], frrext[frwhr].mean())
            frh_mrext[frcur] = np.append(frh_mrext[frcur], frrext[frwhr].max())

            if fr_uflux[frcur] < cflux:
                fr_area[frcur] = len(frwhr[0]) * pix_area
                fr_time[frcur] = dcount
                fr_mlat[frcur] = np.mean(frlat[frwhr[0]])
                fr_mhlcy[frcur] = (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean()
                fr_nhlcy[frcur] = (hlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum()
                fr_sflux[frcur] = (br0[frwhr] * pix_area).sum()
                fr_uflux[frcur] = abs(br0[frwhr] * pix_area).sum()
                fr_rext[frcur] = frrext[frwhr].mean()
                fr_mrext[frcur] = frrext[frwhr].max()

    # Calculate footprint connectivity, and an index of flux rope fieldlines
    #frcmap = np.zeros(frdim, dtype=np.int32)
    #fr_fl = np.zeros(len(fl))
    #for ifl in np.arange(len(hlcy)):
    #    alfpr1 = fl[ifl][0][0]
    #    alfpr2 = fl[ifl][0][-1]
    #    alfpth1 = np.sin(fl[ifl][1][0] + pi/2)
    #    alfpth2 = np.sin(fl[ifl][1][-1] + pi/2)
    #    alfpph1 = fl[ifl][2][0]
    #    alfpph2 = fl[ifl][2][-1]
    #    
    #    alg1 = fl[ifl][0][0] < 1.2
    #    alg2 = fl[ifl][0][-1] < 1.2
    #    
    #    why1 = np.where((frlat_edge < alfpth1) & np.roll((frlat_edge > alfpth1),-1))[0]
    #    whx1 = np.where((frlon_edge < alfpph1) & np.roll((frlon_edge > alfpph1),-1))[0]
    #    why2 = np.where((frlat_edge < alfpth2) & np.roll((frlat_edge > alfpth2),-1))[0]
    #    whx2 = np.where((frlon_edge < alfpph2) & np.roll((frlon_edge > alfpph2),-1))[0]
    #    
    #    ## Map out the pixel connectivity
    #    if (len(whx1 == 1))&(len(why1) == 1)&(frmap[why1,whx1] != 0):
    #        if (len(whx2 == 1))&(len(why2) == 1):
    #            frcmap[why1,whx1] = frmap[why2[0], whx2[0]]
    #            fr_fl[ifl] = frmap[why1, whx1]
    #    if (len(whx2 == 1))&(len(why2) == 1)&(frmap[why2,whx2] != 0):
    #        if (len(whx1 == 1))&(len(why1) == 1):
    #            frcmap[why2,whx2] = frmap[why1[0], whx1[0]]
    #            fr_fl[ifl] = frmap[why2, whx2]

    # Store flux rope and helicity mapping for future use
    frmap0 = np.copy(frmap)
    frhlcy0 = np.copy(frhlcy)

    # Output appropriate arrays to disk
    outfile = netcdf.netcdf_file(outdir + 'fr-' + csfrm + '.nc', 'w')
    outfile.history = 'FRoDO flux rope data'
    outfile.createDimension('lat', frdim[0])
    outfile.createDimension('lon', frdim[1])
    
    out_frmap = outfile.createVariable('frmap', np.int16, ('lat', 'lon'))
    out_frmap[:] = frmap
    out_frmap.units = 'Flux rope footprint label (unitless)'
    
    #out_frcmap = outfile.createVariable('frcmap', np.int16, ('lat', 'lon'))
    #out_frcmap[:] = frcmap
    #out_frcmap.units= 'Flux rope footprint connectivity label (unitless)'
    
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
    if dcount == 0:
        timedel = (time1 - time0)
    else:
        timedel = ((time1 - time0) + timedel) / 2
    timeeta = (nfrm - (dcount+1)) * timedel + time1
    print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta))

    # Advance the timing index
    dcount = dcount + 1

# Create a filter to map those flux ropes that meet the radial extent criteria
fr_hdrext = np.zeros(len(fr_dur), dtype=np.double)
for fri in np.arange(1,len(fr_dur)):
    fr_hdrext[fri] = len(np.where(frh_rext[fri] > 1.5)[0])
fr_rfrg = np.where((fr_hdrext / fr_dur) < 0.5)[0]

# Create a filter to remove flux ropes with only a single day of time history
fr_dfrg = np.where((fr_dur != 1) & (np.isfinite(fr_dur)))[0]

# Merge these into a single filtered index
fr_frg = np.intersect1d(fr_rfrg, fr_dfrg)

# Save this filtered list
outfile = open(outdir + '/hist/fr-frg.pkl', 'wb')
pickle.dump(fr_frg, outfile)
outfile.close()

# Output any completed time-series variables
outfile = open(outdir + '/hist/h-fr-area.pkl', 'wb')
pickle.dump(fr_area, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-time.pkl', 'wb')
pickle.dump(fr_time, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-dur.pkl', 'wb')
pickle.dump(fr_dur, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-mhlcy.pkl', 'wb')
pickle.dump(fr_mhlcy, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-nhlcy.pkl', 'wb')
pickle.dump(fr_nhlcy, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-sflux.pkl', 'wb')
pickle.dump(fr_sflux, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-uflux.pkl', 'wb')
pickle.dump(fr_uflux, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-rext.pkl', 'wb')
pickle.dump(fr_rext, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-mrext.pkl', 'wb')
pickle.dump(fr_mrext, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-mlat.pkl', 'wb')
pickle.dump(fr_mlat, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-fr-tarr.pkl', 'wb')
pickle.dump(tarr, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-area.pkl', 'wb')
pickle.dump(frh_area, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-time.pkl', 'wb')
pickle.dump(frh_time, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-mhlcy.pkl', 'wb')
pickle.dump(frh_mhlcy, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-nhlcy.pkl', 'wb')
pickle.dump(frh_nhlcy, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-sflux.pkl', 'wb')
pickle.dump(frh_sflux, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-uflux.pkl', 'wb')
pickle.dump(frh_uflux, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-rext.pkl', 'wb')
pickle.dump(frh_rext, outfile)
outfile.close()

outfile = open(outdir + '/hist/h-frh-mrext.pkl', 'wb')
pickle.dump(frh_mrext, outfile)
outfile.close()
