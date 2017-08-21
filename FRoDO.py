# FRoDO.py
# Flux Rope Detection and Organization
'''
This set of modules is the Flux Rope Detection and Organization (FRoDO) code.
For ease of use, it can be executed directly from the command line with python3 FRoDO.py
For all additional details, consult the aptly named README file.
'''

# Import libraries
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *

from matplotlib import gridspec
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
import mpl_toolkits.axisartist.floating_axes as floating_axes

import b_sim_netcdf
import os
import glob
import shutil
import re
import datetime
import pickle
import pandas as pd

import numpy as np
from numpy import pi
import scipy
import scipy.stats
from scipy.io import netcdf

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

def prep():
    '''
    Computes magnetic vector potential values in the deVore gauge from input magnetic field data.
    '''

    # A quick announcement
    print('Prepping input data...')

    # Import prep libraries
    from compA import compA

    # Generate a list of files to search through
    files = glob.glob(datdir+bdatprefix+'*.nc')
    files.sort()

    # Cycle through input data
    for file in files:

        # Extract the file frame string
        frmstr = re.split(bdatprefix+'|\.', file)[1]

        # Compute and store magnetic vector potential
        compA.compa(frmstr, datdir, bdatprefix, adatprefix)

def FRoDO():
    '''
    Detects and organizes time histories for magnetic flux ropes.
    '''

    # A quick announcement
    print('Running FRoDO flux rope detection...')

    # Remove any existing files
    if os.path.exists(outdir) : shutil.rmtree(outdir)

    # Create output directories if needed
    if not os.path.exists(outdir) : os.mkdir(outdir)
    if not os.path.exists(outdir + 'hist/') : os.mkdir(outdir + 'hist/')

    # Generate a list of files to search through
    bfrm_list = glob.glob(datdir + bdatprefix + '*.nc')
    bfrm_list.sort()
    afrm_list = glob.glob(datdir + adatprefix + '*.nc')
    afrm_list.sort()
    nfrm = len(bfrm_list)

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
    pix_area = 2. * np.pi * rsun**2 * (np.cos((90-maxlat)*np.pi/180.) - np.cos((90+maxlat)*np.pi/180.)) / (frdim[0] * frdim[1])
    pix_area_pol = 4. * np.pi * rsun**2 / (frdim[0] * frdim[1])

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
    prntend = '\r'

    # Begin cycling through time frames
    for cfrm in bfrm_list:

        # Define some timing
        time0 = datetime.datetime.now()
        csfrm = re.split(bdatprefix+'|\.', cfrm)[1]

        # Read magnetic field data into memory
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
        cjuld = np.double(ctarr.toordinal()+1721425.0)
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
        frlat = np.linspace(-np.sin(maxlat*np.pi/180.)+(1/(2*frdim[0])), np.sin(maxlat*np.pi/180.)-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)
        frlat_pol = np.linspace(-1+(1/(2*frdim[0])), 1-(1/(2*frdim[0])), num=frdim[0], dtype=np.double)

        frlat_edge = frlat - abs(frlat[0] - frlat[1])/2.
        frlon_edge = frlon - abs(frlon[0] - frlon[1])/2.

        # Interpolate the magnetic field array to the output grid
        f0b = scipy.interpolate.interp2d(lons[0,:]*2*pi/360, np.sin(lats[:,0]*2*pi/360), np.rot90(br[:,:,0]), kind='cubic')
        br0 = f0b(frlon, frlat)
        br0_pol = f0b(frlon, frlat_pol)

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
            dtime = 1
        else:
            dtime0 = tarr[dcount] - tarr[dcount - 1]
            dtime = dtime0.days + (dtime0.seconds / 86400.)
        drot = diff_rot(dtime * u.day, np.arcsin(frlat)*180/pi * u.deg, rot_type='snodgrass')
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
                fr_dur = np.append(fr_dur, dtime)
                fr_mhlcy = np.append(fr_mhlcy, (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean())
                fr_nhlcy = np.append(fr_nhlcy, (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum())
                fr_sflux = np.append(fr_sflux, (br0[frwhr] * pix_area).sum())
                fr_uflux = np.append(fr_uflux, abs(br0[frwhr] * pix_area).sum())
                fr_rext = np.append(fr_rext, frrext[frwhr].mean())
                fr_mrext = np.append(fr_mrext, frrext[frwhr].max())

                frh_area.append(np.array([len(frwhr[0]) * pix_area]))
                frh_time.append(np.array([dcount]))
                frh_mhlcy.append(np.array([(frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean()]))
                frh_nhlcy.append(np.array([(frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum()]))
                frh_sflux.append(np.array([(br0[frwhr] * pix_area).sum()]))
                frh_uflux.append(np.array([abs(br0[frwhr] * pix_area).sum()]))
                frh_rext.append(np.array([frrext[frwhr].mean()]))
                frh_mrext.append(np.array([frrext[frwhr].max()]))
            else:
                cflux = abs(br0[frwhr] * pix_area).sum()
                fr_dur[frcur] = fr_dur[frcur] + dtime

                frh_area[frcur] = np.append(frh_area[frcur], len(frwhr[0]) * pix_area)
                frh_time[frcur] = np.append(frh_time[frcur], dcount)
                frh_mhlcy[frcur] = np.append(frh_mhlcy[frcur], (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean())
                frh_nhlcy[frcur] = np.append(frh_nhlcy[frcur], (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum())
                frh_sflux[frcur] = np.append(frh_sflux[frcur], (br0[frwhr] * pix_area).sum())
                frh_uflux[frcur] = np.append(frh_uflux[frcur], abs(br0[frwhr] * pix_area).sum())
                frh_rext[frcur] = np.append(frh_rext[frcur], frrext[frwhr].mean())
                frh_mrext[frcur] = np.append(frh_mrext[frcur], frrext[frwhr].max())

                if fr_uflux[frcur] < cflux:
                    fr_area[frcur] = len(frwhr[0]) * pix_area
                    fr_time[frcur] = dcount
                    fr_mlat[frcur] = np.mean(frlat[frwhr[0]])
                    fr_mhlcy[frcur] = (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).mean()
                    fr_nhlcy[frcur] = (frhlcy[frwhr] * np.abs(br0[frwhr]) * rsun**2 * pix_area).sum()
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
        out_frrext.units = 'Field line maximum radial extent (R_sun)'
        
        out_br0 = outfile.createVariable('br0', np.double, ('lat', 'lon'))
        out_br0[:] = br0
        out_br0.units = 'Radial magnetic flux density at 1.0 R_sun (G)'

        out_br0_pol = outfile.createVariable('br0_pol', np.double, ('lat', 'lon'))
        out_br0_pol[:] = br0_pol
        out_br0_pol.units = 'Radial magnetic flux density at 1.0 R_sun (G)'
        
        out_lat = outfile.createVariable('lat', np.float32, ('lat',))
        out_lat[:] = frlat
        out_lon = outfile.createVariable('lon', np.float32, ('lon',))
        out_lon[:] = frlon

        # Diagnostic readouts!
        time1 = datetime.datetime.now()
        if dcount == 0:
            timedel = (time1 - time0)
        else:
            timedel = ((time1 - time0) + timedel) / 2
        timeeta = (nfrm - (dcount+1)) * timedel + time1
        if dcount == (len(bfrm_list) - 1) : prntend = '\n'
        print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta), end=prntend)

        # Advance the timing index
        dcount = dcount + 1

    # Create a filter to map those flux ropes that meet the radial extent criteria
    dtarr = tarr - np.roll(tarr,1)
    dtarr[0] = datetime.timedelta(1)
    for i in np.arange(len(dtarr)): dtarr[i] = dtarr[i].days
    fr_hdrext = np.zeros(len(fr_dur), dtype=np.double)
    for fri in np.arange(1,len(fr_dur)):
        fr_hdrext[fri] = dtarr[(np.where(frh_rext[fri] > 1.5)[0])].sum()
    fr_rfrg = np.where((fr_hdrext[1:] / fr_dur[1:]) < 0.5)[0] + 1

    # Create a filter to remove flux ropes with less than a single day of time history
    fr_dfrg = np.where((fr_dur[1:] > 1) & (np.isfinite(fr_dur[1:])))[0]+1

    # Merge these into a single filtered index
    fr_frg = np.intersect1d(fr_rfrg, fr_dfrg)

    # Save this filtered list
    outfile = open(outdir + '/hist/fr-frg.pkl', 'wb')
    pickle.dump(fr_frg, outfile)
    outfile.close()

    # Output fieldline helicity threshold values
    outfile = open(outdir + '/hist/ethreshs.pkl', 'wb')
    pickle.dump(ethreshs, outfile)
    outfile.close()

    outfile = open(outdir + '/hist/sthreshs.pkl', 'wb')
    pickle.dump(sthreshs, outfile)
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

def erupt():
    '''
    Detects and flags erupting flux rope signatures.
    '''

    # A quick announcement
    print('Running FRoDO flux rope eruption detection...')

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
    prntend = '\r'

    for cfrm in bfrm_list:

        # Define some timing
        time0 = datetime.datetime.now()
        csfrm = re.split(bdatprefix+'|\.', cfrm)[1]

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
        cjuld = np.double(ctarr.toordinal()+1721425.0)
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
                    if len(fpwhr[ifpfr][1]) > fpmaxarea:
                        fpmaxarea = len(fpwhr[ifpfr][1])
                        fpmaxarea_reg = int(fpreg[ifpfr])
            if (fpmaxarea_reg != 0):
                fr_etarr = np.append(fr_etarr, dcount)
                fr_efpt = np.append(fr_efpt, fpmaxarea_reg)

        # Diagnostic readouts!
        time1 = datetime.datetime.now()
        if dcount == 0:
            timedel = (time1 - time0)
        else:
            timedel = ((time1 - time0) + timedel) / 2
        timeeta = (nfrm - (dcount+1)) * timedel + time1
        if dcount == (len(bfrm_list) - 1) : prntend = '\n'
        print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta), end=prntend)

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

def plot():
    '''
    Creates a standard set of plot outputs for detected flux ropes.
    '''

    # A quick announcement
    print('Running plotting routines...')

    # Define color tables
    import palettable
    cols = (palettable.colorbrewer.get_map('Paired', 'Qualitative', 12)).mpl_colors

    # Define some plotting parameters
    fscale = 1.0            # Relative figure scale
    fnapp = ''              # Label to append to plot filename
    defcol = 'k'            # Default primary color
    defcol2 = cols[5]       # Default secondary color
    erptcol = '#ff7f0e'     # Erupting plotting color
    erptcmap = 'Oranges'    # Erupting color map
    nerptcol = '#1f77B4'    # Non-erupting plotting color
    nerptcmap = 'Blues'     # Non-erupting color map
    legfsz = 8              # Legend font size
    gsize = 1e20            # Butterfly glyph scaling size

    # Remove any existing files
    if os.path.exists('plt') : shutil.rmtree('plt')

    # Create output directories if needed
    if not os.path.exists('plt') : os.mkdir('plt')

    # Define a quick pre-built function for radial plotting
    # Reverse-engineered from a matplotlib example script
    def setup_polaraxis(fig, rect, radlim):
        """
        Sometimes, things like axis_direction need to be adjusted.
        """

        # rotate a bit for better orientation
        tr_rotate = Affine2D().translate(90, 0)

        # scale degree to radians
        tr_scale = Affine2D().scale(np.pi/180., 1.)

        tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

        grid_locator1 = angle_helper.LocatorHMS(8)
        tick_formatter1 = angle_helper.FormatterHMS()

        grid_locator2 = MaxNLocator(3)

        ra0, ra1 = -90, 90
        cz0, cz1 = 0, radlim
        grid_helper = floating_axes.GridHelperCurveLinear(
            tr, extremes=(ra0, ra1, cz0, cz1),
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
            tick_formatter1=None,
            tick_formatter2=None)

        ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
        fig.add_subplot(ax1)

        # adjust axis
        ax1.axis["left"].set_axis_direction("bottom")
        ax1.axis["right"].set_axis_direction("top")

        ax1.axis["bottom"].set_visible(False)
        ax1.axis["top"].set_axis_direction("bottom")
        ax1.axis["top"].toggle(ticklabels=True, label=True)
        ax1.axis["top"].major_ticklabels.set_axis_direction("top")
        ax1.axis["top"].label.set_axis_direction("top")

        ax1.axis["left"].label.set_text('Number [bin$^{-1}$]')
        ax1.axis["top"].label.set_text(r'$\theta$ [degrees]')

        # create a parasite axes whose transData in RA, cz
        aux_ax = ax1.get_aux_axes(tr)

        aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
        ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
        # drawn twice, and possibly over some other
        # artists. So, we decrease the zorder a bit to
        # prevent this.

        return ax1, aux_ax

    # Read time histories
    infile = open(outdir + '/hist/h-fr-area.pkl', 'rb')
    fr_area = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-time.pkl', 'rb')
    fr_time = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-dur.pkl', 'rb')
    fr_dur = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mhlcy.pkl', 'rb')
    fr_mhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-nhlcy.pkl', 'rb')
    fr_nhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-sflux.pkl', 'rb')
    fr_sflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-uflux.pkl', 'rb')
    fr_uflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-rext.pkl', 'rb')
    fr_rext = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mrext.pkl', 'rb')
    fr_mrext = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mlat.pkl', 'rb')
    fr_mlat = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-tarr.pkl', 'rb')
    tarr = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-area.pkl', 'rb')
    frh_area = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-time.pkl', 'rb')
    frh_time = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-mhlcy.pkl', 'rb')
    frh_mhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-nhlcy.pkl', 'rb')
    frh_nhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-sflux.pkl', 'rb')
    frh_sflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-uflux.pkl', 'rb')
    frh_uflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-rext.pkl', 'rb')
    frh_rext = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-frh-mrext.pkl', 'rb')
    frh_mrext = pickle.load(infile)
    infile.close()

    # Read eruption labels
    infile = open(outdir + '/hist/fr-efpt.pkl', 'rb')
    fr_efpt = pickle.load(infile)
    infile.close()
    fr_efpt = fr_efpt.astype(np.int)

    infile = open(outdir + '/hist/fr-elab.pkl', 'rb')
    fr_elab = pickle.load(infile)
    infile.close()
    fr_elab = fr_elab.astype(np.int)

    infile = open(outdir + '/hist/fr-etarr.pkl', 'rb')
    fr_etarr = pickle.load(infile)
    infile.close()

    # Read radial extent and duration filtered index
    infile = open(outdir + '/hist/fr-frg.pkl', 'rb')
    fr_frg = pickle.load(infile)
    infile.close()
    fr_frg = fr_frg.astype(np.int)

    ## Create an index of non-erupting flux ropes
    regmask = np.ones(len(fr_area), np.bool)
    regmask[fr_elab] = 0
    fr_nelab = np.where(regmask)[0]
    del regmask

    # Merge the set of erupting structures with the list of confirmed flux ropes from radial extent
    fr_elab = fr_elab[np.in1d(fr_elab, fr_frg)]
    fr_nelab = fr_nelab[np.in1d(fr_nelab, fr_frg)]
    fr_etarr = (fr_etarr[np.in1d(fr_efpt, fr_frg)]).astype(np.int)
    fr_efpt = fr_efpt[np.in1d(fr_efpt, fr_frg)]

    # Define time arrays
    netarr = []
    etarr = []
    erptarr = []
    nearr = (fr_time[fr_nelab]).astype(np.int)
    earr = (fr_time[fr_elab]).astype(np.int)
    for i in np.arange(len(nearr)):
        netarr.append(tarr[nearr[i]])
    for i in np.arange(len(earr)):
        etarr.append(tarr[earr[i]])
    for i in np.arange(len(fr_etarr)):
        erptarr.append(tarr[fr_etarr[i]])

    # Compute the flux and helicity ejection rates
    tarr = np.array(tarr)
    ejr_uflux = np.zeros(len(tarr), dtype=np.double)
    ejr_nhlcy = np.zeros(len(tarr), dtype=np.double)
    ejr_nerupt = np.zeros(len(tarr), dtype=np.double)
    for i in fr_efpt:
        ejr_uflux[fr_time[i].astype(np.int)] = ejr_uflux[fr_time[i].astype(np.int)] + fr_uflux[i]
        ejr_nhlcy[fr_time[i].astype(np.int)] = ejr_nhlcy[fr_time[i].astype(np.int)] + abs(fr_nhlcy[i])
        ejr_nerupt[fr_time[i].astype(np.int)] = ejr_nerupt[fr_time[i].astype(np.int)] + 1

    dtarr = pd.to_datetime(np.array(pd.date_range(tarr[0].date(), tarr[-1].date(), freq='1D')))
    dejr_uflux = np.zeros(len(dtarr), dtype=np.double)
    dejr_nhlcy = np.zeros(len(dtarr), dtype=np.double)
    dejr_nerupt = np.zeros(len(dtarr), dtype=np.double)
    for i in np.arange(len(dtarr)):
        if i != len(dtarr)-1:
           dejr_uflux[i] = ejr_uflux[(tarr >= dtarr[i]) & (tarr < dtarr[i+1])].sum()
           dejr_nhlcy[i] = ejr_nhlcy[(tarr >= dtarr[i]) & (tarr < dtarr[i+1])].sum()
           dejr_nerupt[i] = ejr_nerupt[(tarr >= dtarr[i]) & (tarr < dtarr[i+1])].sum()
        else:
           dejr_uflux[i] = ejr_uflux[tarr >= dtarr[i]].sum()
           dejr_nhlcy[i] = ejr_nhlcy[tarr >= dtarr[i]].sum()
           dejr_nerupt[i] = ejr_nerupt[tarr >= dtarr[i]].sum()

    # Generate a persistence array for flux ropes
    ofiles = glob.glob(outdir + 'fr-*.nc')
    ofiles.sort()
    fr_pmap = np.zeros(frdim, dtype=np.double)
    fcount = 0
    for file in ofiles:
        if fcount != 0:
            dtime0 = tarr[fcount] - tarr[fcount - 1]
            dtime = dtime0.days + (dtime0.seconds / 86400.)
        else:
            dtime = 1.
        infile = netcdf.netcdf_file(file, 'r')
        fr_pmap = fr_pmap + (infile.variables['frmap'][:,:] != 0) * dtime
        infile.close()
        fcount = fcount + 1

    # Begin plotting

    # Radial graph of eruption mean latitudes
    fpmlatne = np.histogram(np.arcsin(fr_mlat[fr_nelab])*180./np.pi, bins=90, range=[-90,90], density=False)
    fpmlate = np.histogram(np.arcsin(fr_mlat[fr_elab])*180./np.pi, bins=90, range=[-90,90], density=False)
    whg = (180.) / len(fpmlatne[0])
    rscl = np.array([fpmlatne[0].max(), fpmlate[0].max()]).max() * 1.2

    f = figure(figsize=fscale*np.array([3.35,2.45]))
    ax, aax = setup_polaraxis(f, 111, rscl)
    aax.plot(fpmlatne[1][0:-1], fpmlatne[0], label='Non-erupting', color=nerptcol)
    aax.plot(fpmlate[1][0:-1], fpmlate[0], label='Erupting', color=erptcol)
    aax.legend(loc=(0.65,1.00))
    ax.grid()
    tight_layout()
    f.savefig('plt/fr-mlats'+fnapp+'.pdf')

    # Butterfly digram of flux rope footprints
    f = figure(figsize=fscale*np.array([7,4]))
    gs = gridspec.GridSpec(2, 1,
                   width_ratios=[1],
                   height_ratios=[1,1]
                   )
    ax2 = f.add_subplot(gs[0])
    ax1 = f.add_subplot(gs[1])
    sctne = ax1.scatter(netarr, np.arcsin(fr_mlat[fr_nelab])*180./np.pi, c=fr_nhlcy[fr_nelab], cmap='RdBu_r', s=fr_area[fr_nelab]/gsize,edgecolors='None',vmin=-0.5e43,vmax=0.5e43, alpha=1.0)
    scte = ax2.scatter(etarr, np.arcsin(fr_mlat[fr_elab])*180./np.pi, c=fr_nhlcy[fr_elab], cmap='RdBu_r', s=fr_area[fr_elab]/gsize,edgecolors='None',vmin=-0.5e43,vmax=0.5e43, alpha=1.0)
    cb = colorbar(sctne,ax=ax1, extend='both', label='Non-erupting H [Mx$^2$]', fraction=0.05)
    cb2 = colorbar(scte,ax=ax2, extend='both', label='Erupting H [Mx$^2$]', fraction=0.05)
    ax1.set_ylim([-90,90])
    ax1.set_axisbelow(True)
    ax1.set_xlim([tarr[0], tarr[-1]])
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylim([-90,90])
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Latitude')
    ax2.set_ylabel('Latitude')
    ax2.set_axisbelow(True)
    ax2.set_xlim([tarr[0], tarr[-1]])
    ax1.grid()
    ax2.grid()
    f.autofmt_xdate() # Optional toggle to sort out overlapping dates
    tight_layout()
    f.savefig('plt/fr-bfly-nshlcy'+fnapp+'.pdf')

    # Butterfly digram of erupting flux ropes
    f, (ax1) = subplots(figsize=fscale*np.array([7,2]))
    sctne = ax1.scatter(erptarr, np.arcsin(fr_mlat[fr_efpt])*180./np.pi, c=fr_nhlcy[fr_efpt], cmap='RdBu_r', s=fr_area[fr_efpt]/gsize,edgecolors='None',vmin=-0.5e43,vmax=0.5e43, alpha=1.0)
    cb = colorbar(sctne,ax=ax1, extend='both', label='Erupting H [Mx$^2$]', fraction=0.05)
    ax1.set_ylim([-90,90])
    ax1.set_axisbelow(True)
    ax1.set_xlim([tarr[0], tarr[-1]])
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Latitude')
    ax1.grid()
    #f.autofmt_xdate() # Optional toggle to sort out overlapping dates
    tight_layout()
    f.savefig('plt/fr-bfly-nshlcy-erupt'+fnapp+'.pdf')

    # A few hexbin scatter plots

    # Scatter of duration and net helicity magnitude
    if len(fr_nelab) != 0:
        ne_coef = np.polyfit(np.log10(fr_dur[fr_nelab]), np.log10(abs(fr_nhlcy[fr_nelab])),1)
    else: ne_coef = np.array([0,0])
    ne_polyn = np.poly1d(ne_coef)
    ne_xs = np.linspace(0.1,2.8)
    ne_ys = ne_polyn(ne_xs)

    if len(fr_nelab) != 0:
        e_coef = np.polyfit(np.log10(fr_dur[fr_elab]), np.log10(abs(fr_nhlcy[fr_elab])),1)
    else: ecoef = np.array([0,0])
    e_polyn = np.poly1d(e_coef)
    e_xs = np.linspace(0.1,2.8)
    e_ys = e_polyn(e_xs)

    f, (ax1,ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
    hb1 = ax1.hexbin(fr_dur[fr_nelab], abs(fr_nhlcy[fr_nelab]), xscale='log', yscale='log', gridsize=25, extent=[0,3,40,44], cmap=nerptcmap, linewidths=0)
    hb2 = ax2.hexbin(fr_dur[fr_elab], abs(fr_nhlcy[fr_elab]), xscale='log', yscale='log', gridsize=25, extent=[0,3,40,44], cmap=erptcmap, linewidths=0)
    ax2.set_xlabel('Duration [days]')
    ax2.set_ylabel('$|$H$|$ [Mx$^2$]')
    ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+'} x^{'+'%.03f'%ne_coef[0]+'}$')
    ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+'} x^{'+'%.03f'%e_coef[0]+'}$')
    if len(fr_nelab) != 0 : leg1 = ax1.legend(loc=4, fontsize=legfsz)
    if len(fr_elab) != 0 : leg2 = ax2.legend(loc=4, fontsize=legfsz)
    ax1.set_xlim(1e0,1e3)
    ax1.set_ylim(1e40,1e44)
    ax2.set_xlim(1e0,1e3)
    ax2.set_ylim(1e40,1e44)
    ax1.xaxis.set_ticklabels([])
    if len(fr_nelab) == 0 : hb1.set_clim(0,1)
    if len(fr_elab) == 0 : hb2.set_clim(0,1)
    colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
    colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
    tight_layout()
    savefig('plt/fr-sct-dur-nhlcy'+fnapp+'.pdf')

    # Scatter of duration and unsigned magnetic flux
    if len(fr_nelab) != 0:
        ne_coef = np.polyfit(np.log10(fr_dur[fr_nelab]), np.log10(fr_uflux[fr_nelab]),1)
    else: ne_coef = np.array([0,0])
    ne_polyn = np.poly1d(ne_coef)
    ne_xs = np.linspace(0.1,2.8)
    ne_ys = ne_polyn(ne_xs)

    if len(fr_elab) != 0:
        e_coef = np.polyfit(np.log10(fr_dur[fr_elab]), np.log10(fr_uflux[fr_elab]),1)
    else: e_coef = np.array([0,0])
    e_polyn = np.poly1d(e_coef)
    e_xs = np.linspace(0.1,2.8)
    e_ys = e_polyn(e_xs)

    f, (ax1,ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
    hb1 = ax1.hexbin(fr_dur[fr_nelab], fr_uflux[fr_nelab], xscale='log', yscale='log', gridsize=25, extent=[0,3,18,23], cmap=nerptcmap, linewidths=0)
    hb2 = ax2.hexbin(fr_dur[fr_elab], fr_uflux[fr_elab], xscale='log', yscale='log', gridsize=25, extent=[0,3,18,23], cmap=erptcmap, linewidths=0)
    ax2.set_xlabel('Duration [days]')
    ax2.set_ylabel('$|\Phi_m|$ [Mx]')
    ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+  '} x^{'+'%.03f'%ne_coef[0]+'}$')
    ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+  '} x^{'+'%.03f'%e_coef[0]+'}$')
    if len(fr_nelab) != 0 : leg1 = ax1.legend(loc=4, fontsize=legfsz)
    if len(fr_elab) != 0: leg2 = ax2.legend(loc=4, fontsize=legfsz)
    ax1.set_xlim(1e0,1e3)
    ax1.set_ylim(1e18,1e23)
    ax2.set_xlim(1e0,1e3)
    ax2.set_ylim(1e18,1e23)
    ax1.xaxis.set_ticklabels([])
    if len(fr_nelab) == 0 : hb1.set_clim(0,1)
    if len(fr_elab) == 0 : hb2.set_clim(0,1)
    colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
    colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
    tight_layout()
    savefig('plt/fr-sct-dur-uflux'+fnapp+'.pdf')

    # Scatter of unsigned magnetic flux and net helicity magnitude
    if len(fr_nelab) != 0:
        ne_coef = np.polyfit(np.log10(fr_uflux[fr_nelab]), np.log10(abs(fr_nhlcy[fr_nelab])),1)
    else: ne_coef = np.array([0,0])
    ne_polyn = np.poly1d(ne_coef)
    ne_xs = np.linspace(18.1,22.8)
    ne_ys = ne_polyn(ne_xs)

    if len(fr_elab) != 0:
        e_coef = np.polyfit(np.log10(fr_uflux[fr_elab]), np.log10(abs(fr_nhlcy[fr_elab])),1)
    else: e_coef = np.array([0,0])
    e_polyn = np.poly1d(e_coef)
    e_xs = np.linspace(18.1,22.8)
    e_ys = e_polyn(e_xs)

    f, (ax1, ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
    hb1=ax1.hexbin(fr_uflux[fr_nelab], abs(fr_nhlcy[fr_nelab]), xscale='log', yscale='log', gridsize=25, cmap=nerptcmap, extent=[18,23,40,44], linewidths=0)
    hb2=ax2.hexbin(fr_uflux[fr_elab], abs(fr_nhlcy[fr_elab]), xscale='log', yscale='log', gridsize=25, cmap=erptcmap, extent=[18,23,40,44], alpha=1.0, linewidths=0)
    #ax1.grid()
    #ax2.grid()
    ax2.set_xlabel(r'$|\Phi_m|$ [Mx]')
    ax2.set_ylabel(r'$|$H$|$ [Mx$^2$]')
    ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+  '} x^{'+'%.03f'%ne_coef[0]+'}$')
    ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+  '} x^{'+'%.03f'%e_coef[0]+'}$')
    if len(fr_nelab) != 0 : leg1 = ax1.legend(loc=4, fontsize=legfsz)
    if len(fr_elab) != 0: leg2 = ax2.legend(loc=4, fontsize=legfsz)
    ax1.set_xlim(1e18,1e23)
    ax1.set_ylim(1e40,1e44)
    ax2.set_xlim(1e18,1e23)
    ax2.set_ylim(1e40,1e44)
    ax1.xaxis.set_ticklabels([])
    if len(fr_nelab) == 0 : hb1.set_clim(0,1)
    if len(fr_elab) == 0 : hb2.set_clim(0,1)
    colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
    colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
    tight_layout()
    f.savefig('plt/fr-sct-uflux-nhlcy'+fnapp+'.pdf')

    # Helicity and unsigned flux ejection rates
    f = figure(figsize=fscale*np.array([3.35,5]))
    gs = gridspec.GridSpec(3, 1,
                   width_ratios=[1],
                   height_ratios=[1,2,2]
                   )
    ax0 = f.add_subplot(gs[0])
    ax1 = f.add_subplot(gs[1])
    ax2 = f.add_subplot(gs[2])

    ax0.plot(tarr, ejr_nerupt)
    ax0.plot(dtarr, scipy.convolve(dejr_nerupt, np.ones(7)/7., mode='same'), label='7-day', color=defcol2)
    ax0.set_ylabel('$N_{erupt}$')

    ax1.plot(tarr, ejr_uflux, label='Original', color=defcol)
    ax1.plot(dtarr, scipy.convolve(dejr_uflux, np.ones(7)/7., mode='same'), label='7-day', color=defcol2)
    ax2.plot(tarr, ejr_nhlcy, label='Original', color=defcol)
    ax2.plot(dtarr, scipy.convolve(dejr_nhlcy, np.ones(7)/7., mode='same'), label='7-day', color=defcol2)

    ax2.set_xlabel('Date')
    ax1.set_ylabel(r'$\Phi_m$ [Mx]')
    ax2.set_ylabel(r'H [Mx$^2$]')
    ax0.xaxis.set_ticklabels([])
    ax1.xaxis.set_ticklabels([])
    for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    ax1.legend()
    #f.autofmt_xdate() # Optional toggle to sort out overlapping dates
    tight_layout()
    savefig('plt/fr_ejr'+fnapp+'.pdf')

    # Flux rope persistence map
    f, (ax) = subplots(1, figsize=fscale*np.array([3.35,2.0]))
    im = ax.imshow(fr_pmap, extent=[0,360,-1,1], cmap='Greys', aspect='auto')
    colorbar(im, label='Flux rope persistence [days]')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Sine latitude')
    tight_layout()
    savefig('plt/fr_pmap'+fnapp+'.pdf')

def stats(tex=False):
    '''
    Computes a series of statistics for erupting and non-erupting magnetic flux ropes.
    '''

    # A quick announcement
    print('Calculating FRoDO flux rope statistics...')

    # Remove any existing files
    if os.path.exists(outdir + 'stats.txt') : os.remove(outdir + 'stats.txt')
    if os.path.exists(outdir + 'stats-tex.tex') : os.remove(outdir + 'stats-tex.tex')

    # Read data
    infile = open(outdir + '/hist/h-fr-area.pkl', 'rb')
    fr_area = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-time.pkl', 'rb')
    fr_time = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-dur.pkl', 'rb')
    fr_dur = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mhlcy.pkl', 'rb')
    fr_mhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-nhlcy.pkl', 'rb')
    fr_nhlcy = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-sflux.pkl', 'rb')
    fr_sflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-uflux.pkl', 'rb')
    fr_uflux = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-rext.pkl', 'rb')
    fr_rext = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mrext.pkl', 'rb')
    fr_mrext = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-mlat.pkl', 'rb')
    fr_mlat = pickle.load(infile)
    infile.close()

    infile = open(outdir + '/hist/h-fr-tarr.pkl', 'rb')
    tarr = pickle.load(infile)
    infile.close()

    # Read eruption labels
    infile = open(outdir + '/hist/fr-efpt.pkl', 'rb')
    fr_efpt = pickle.load(infile)
    infile.close()
    fr_efpt = fr_efpt.astype(np.int)

    infile = open(outdir + '/hist/fr-elab.pkl', 'rb')
    fr_elab = pickle.load(infile)
    infile.close()
    fr_elab = fr_elab.astype(np.int)

    infile = open(outdir + '/hist/fr-etarr.pkl', 'rb')
    fr_etarr = pickle.load(infile)
    infile.close()

    # Read radial extent and duration filtered index
    infile = open(outdir + '/hist/fr-frg.pkl', 'rb')
    fr_frg = pickle.load(infile)
    infile.close()
    fr_frg = fr_frg.astype(np.int)

    ## Create an index of non-erupting flux ropes
    regmask = np.ones(len(fr_area), np.bool)
    regmask[fr_elab] = 0
    fr_nelab = np.where(regmask)[0]
    del regmask

    # Merge the set of erupting structures with the list of confirmed flux ropes from radial extent
    fr_elab = fr_elab[np.in1d(fr_elab, fr_frg)]
    fr_nelab = fr_nelab[np.in1d(fr_nelab, fr_frg)]
    fr_efpt = fr_efpt[np.in1d(fr_efpt, fr_frg)]

    # Define time arrays
    netarr = []
    etarr = []
    nearr = (fr_time[fr_nelab]).astype(np.int)
    earr = (fr_time[fr_elab]).astype(np.int)
    for i in np.arange(len(nearr)):
        netarr.append(tarr[nearr[i]])
    for i in np.arange(len(earr)):
        etarr.append(tarr[earr[i]])

    # Open file for writing
    f = open(outdir + 'stats.txt', 'w')

    # Calculate and print some statistics

    f.write('Mean erupting unsigned net helicity (Mx^2) ' + '%1.2E'%abs(fr_nhlcy[fr_elab]).mean() + '\n')
    f.write('Mean non-erupting unsigned net helicity (Mx^2) ' + '%1.2E'%abs(fr_nhlcy[fr_nelab]).mean() + '\n')
    f.write('Standard deviation of erupting unsigned net helicity (Mx^2) ' + '%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_elab])) + '\n')
    f.write('Standard deviation non-erupting unsigned net helicity (Mx^2) ' + '%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_nelab])) + '\n')

    f.write('Mean erupting unsigned magnetic flux (Mx) ' + '%1.2E'%abs(fr_uflux[fr_elab]).mean() + '\n')
    f.write('Mean non-erupting unsigned magnetic flux (Mx) ' + '%1.2E'%abs(fr_uflux[fr_nelab]).mean() + '\n')
    f.write('Standard deviation erupting unsigned magnetic flux (Mx) ' + '%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_elab])) + '\n')
    f.write('Standard deviation non-erupting unsigned magnetic flux (Mx) ' + '%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_nelab])) + '\n')

    f.write('Mean erupting footprint area (cm^2) ' + '%1.2E'%abs(fr_area[fr_elab]).mean() + '\n')
    f.write('Mean non-erupting footprint area (cm^2) ' + '%1.2E'%abs(fr_area[fr_nelab]).mean() + '\n')
    f.write('Standard deviation erupting footprint area (cm^2) ' + '%1.2E'%scipy.stats.tstd(abs(fr_area[fr_elab])) + '\n')
    f.write('Standard deviation non-erupting footprint area (cm^2) ' + '%1.2E'%scipy.stats.tstd(abs(fr_area[fr_nelab])) + '\n')

    f.write('Mean erupting duration (days) ' + '%2.1f'%abs(fr_dur[fr_elab]).mean() + '\n')
    f.write('Mean non-erupting duration (days) ' + '%2.1f'%abs(fr_dur[fr_nelab]).mean() + '\n')
    f.write('Standard deviation erupting duration (days) ' + '%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_elab])) + '\n')
    f.write('Standard deviation non-erupting duration (days) ' + '%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_nelab])) + '\n')

    # Calculate the t-test for the erupting and non-erupting sets
    # The independent t-test should be used here, as these are separate sets
    # Note that Welchs t-test is used here, as we have differing sample sizes, with unsure equality of population variance.
    # This calculation returns the calculated t-statistic value, and the two-tailed p-value.

    f.write('T-test t-statistics and two-tailed p-values:' + '\n')

    t1 = scipy.stats.ttest_ind(abs(fr_nhlcy[fr_elab]), abs(fr_nhlcy[fr_nelab]),equal_var=False)
    f.write('T-test for unsigned net helicity: ' + '%2.1f'%t1.statistic + ' ' + '%1.2E'%t1.pvalue + '\n')

    t2 = scipy.stats.ttest_ind(fr_uflux[fr_elab], fr_uflux[fr_nelab],equal_var=False)
    f.write('T-test for unsigned magnetic flux: ' + '%2.1f'%t2.statistic + ' ' + '%1.2E'%t2.pvalue + '\n')

    t3 = scipy.stats.ttest_ind(fr_area[fr_elab], fr_area[fr_nelab],equal_var=False)
    f.write('T-test for footprint area: ' + '%2.1f'%t3.statistic + ' ' + '%1.2E'%t3.pvalue + '\n')

    t4 = scipy.stats.ttest_ind(fr_dur[fr_elab], fr_dur[fr_nelab],equal_var=False)
    f.write('T-test for duration: ' + '%2.2f'%t4.statistic + ' ' + '%2.1E'%t4.pvalue + '\n')

    # Calculate linear fits for scatter distributions, along with appropriate statistics.

    f.write('Spearman correlation and p-values:' + '\n')

    # Duration and unsigned net helicity
    s1n = scipy.stats.spearmanr(fr_dur[fr_nelab], abs(fr_nhlcy[fr_nelab]))
    f.write('Spearman rank-order calculation for non-erupting duration and unsigned net helicity: ' + '%1.2f'%s1n.correlation + ' ' +  '%1.1E'%s1n.pvalue + '\n')

    s1e = scipy.stats.spearmanr(fr_dur[fr_elab], abs(fr_nhlcy[fr_elab]))
    f.write('Spearman rank-order calculation for erupting duration and unsigned net helicity: ' + '%1.2f'%s1e.correlation + ' ' + '%1.1E'%s1e.pvalue + '\n')

    # Duration and unsigned magnetic flux
    s2n = scipy.stats.spearmanr(fr_dur[fr_nelab], fr_uflux[fr_nelab])
    f.write('Spearman rank-order calculation for non-erupting duration and unsigned magnetic flux: ' + '%1.2f'%s2n.correlation + ' ' + '%1.1E'%s2n.pvalue + '\n')

    s2e = scipy.stats.spearmanr(fr_dur[fr_elab], fr_uflux[fr_elab])
    f.write('Spearman rank-order calculation for erupting duration and unsigned magnetic flux: ' + '%1.2f'%s2e.correlation + ' ' + '%1.1E'%s2e.pvalue + '\n')

    # Unsigned magnetic flux and unsigned net helicity
    s3n = scipy.stats.spearmanr(fr_uflux[fr_nelab], abs(fr_nhlcy[fr_nelab]))
    f.write('Spearman rank-order calculation for non-erupting unsigned magnetic flux  and unsigned net helicity: ' + '%1.2f'%s3n.correlation + ' ' + '%1.1E'%s3n.pvalue + '\n')

    s3e = scipy.stats.spearmanr(fr_uflux[fr_elab], abs(fr_nhlcy[fr_elab]))
    f.write('Spearman rank-order calculation for erupting unsigned magnetic flux  and unsigned net helicity: ' + '%1.2f'%s3e.correlation + ' ' + '%1.1E'%s3e.pvalue + '\n')

    # All good things...
    f.close()

    # For an optional flag on execution, output the same statistics, but in a friendly LaTeX table format...

    if tex:
        # Open file for writing
        f = open(outdir + 'stats-tex.tex', 'w')

        # Write table of flux rope fit statistics
        f.write(r'\begin{deluxetable}{cccll}'+'\n')
        f.write(r'\tablecaption{Flux rope statistics}'+'\n')
        f.write(r'\tablecolumns{5}'+'\n')
        f.write(r'\tablenum{1}'+'\n')
        f.write(r'\label{tab:fluxrope_stats}'+'\n')
        f.write(r'\tablewidth{0pt}'+'\n')
        f.write(r'\tablehead{'+'\n')
        f.write(r'\colhead{} & \multicolumn{2}{c}{Quantity} & \multicolumn{2}{c}{Spearman}\\'+'\n')
        f.write(r'\cline{2-3} \cline{4-5}\\'+'\n')
        f.write(r'\colhead{E/NE} & \colhead{1} & \colhead{2} & \colhead{cc} & \colhead{$p$-value}}'+'\n')
        f.write(r'\startdata'+'\n')
        f.write(r'\decimals'+'\n')

        # Duration and unsigned net helicity
        f.write('NE & t & $|$H$|$ & ' + '%1.2f'%s1n.correlation + ' & ' + ('%1.1E'%s1n.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s1n.pvalue)[-3:] + '}$ \\\ \n')
        f.write('E & t & $|$H$|$ & ' + '%1.2f'%s1e.correlation + ' & ' + ('%1.1E'%s1e.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s1e.pvalue)[-3:] + '}$ \\\ \n')

        # Duration and unsigned magnetic flux
        f.write('NE & t & $|\Phi_m|$ & ' + '%1.2f'%s2n.correlation + ' & ' + ('%1.1E'%s2n.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s2n.pvalue)[-3:] + '}$ \\\ \n')
        f.write('E & t & $|\Phi_m|$ & ' + '%1.2f'%s2e.correlation + ' & ' + ('%1.1E'%s2e.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s2e.pvalue)[-3:] + '}$ \\\ \n')

        # Unsigned magnetic flux and unsigned net helicity
        f.write('NE & $|\Phi_m|$ & $|$H$|$ & ' + '%1.2f'%s3n.correlation + ' & ' + ('%1.1E'%s3n.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s3n.pvalue)[-3:] + '}$ \\\ \n')
        f.write('E & $|\Phi_m|$ & $|$H$|$ & ' + '%1.2f'%s3e.correlation + ' & ' + ('%1.1E'%s3e.pvalue)[0:3] + '$\\times 10^{' + ('%1.1E'%s3e.pvalue)[-3:] + '}$ \\\ \n')

        f.write('\enddata' + '\n')
        f.write('\end{deluxetable}' + '\n')

        f.write('\n')
        f.write('%%%' + '\n')
        f.write('\n')

        # Write table of flux rope statistics

        f.write(r'\begin{deluxetable*}{rllll}' + '\n')
        f.write(r'\tablecaption{Mean flux rope parameters}' + '\n')
        f.write(r'\tablecolumns{5}' + '\n')
        f.write(r'\tablenum{2}' + '\n')
        f.write(r'\label{tab:fluxrope_params}' + '\n')
        f.write(r'\tablewidth{0pt}' + '\n')
        f.write(r'\tablehead{' + '\n')
        f.write(r'\colhead{Quantity} & \colhead{Erupting} & \colhead{Non-erupting} & \colhead{t-statistic} & \colhead{$p$-value}}' + '\n')
        f.write(r'\decimals' + '\n')
        f.write(r'\startdata' + '\n')

        # Net helicity magnitude
        f.write('$|\\textrm{H}|$ (Mx$^2$) & ' + ('%1.2E'%abs(fr_nhlcy[fr_elab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_nhlcy[fr_elab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_elab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_elab])))[-2:] + '}$ & ' + ('%1.2E'%abs(fr_nhlcy[fr_nelab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_nhlcy[fr_nelab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_nelab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_nelab])))[-2:] + '}$ & ' + '%2.1f'%t1.statistic + ' & ' + ('%1.2E'%t1.pvalue)[0:4] + ' $\\times 10^{' + ('%1.2E'%t1.pvalue)[-3:] + '}$' + '\\\ \n')

        # Unsigned magnetic flux
        f.write('${\Phi}_m$ (Mx) & ' + ('%1.2E'%abs(fr_uflux[fr_elab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_uflux[fr_elab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_elab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_elab])))[-2:] + '}$ & ' + ('%1.2E'%abs(fr_uflux[fr_nelab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_uflux[fr_nelab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_nelab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_nelab])))[-2:] + '}$ & ' + '%2.1f'%t2.statistic + ' & ' + ('%1.2E'%t2.pvalue)[0:4] + ' $\\times 10^{' + ('%1.2E'%t2.pvalue)[-3:] + '}$' + '\\\ \n')

        # Footprint area
        f.write('${\\textrm{A}}$ (cm$^2$) & ' + ('%1.2E'%abs(fr_area[fr_elab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_area[fr_elab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_area[fr_elab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_area[fr_elab])))[-2:] + '}$ & ' + ('%1.2E'%abs(fr_area[fr_nelab]).mean())[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%abs(fr_area[fr_nelab]).mean())[-2:] + '}$ $\pm$ ' + ('%1.2E'%scipy.stats.tstd(abs(fr_area[fr_nelab])))[0:4] + ' $ ' + '\\times 10^{' + ('%1.2E'%scipy.stats.tstd(abs(fr_area[fr_nelab])))[-2:] + '}$ & ' + '%2.1f'%t3.statistic + ' & ' + ('%1.2E'%t3.pvalue)[0:4] + ' $\\times 10^{' + ('%1.2E'%t3.pvalue)[-3:] + '}$' + '\\\ \n')

        # Duration
        f.write('${\\tau}$ (days) & ' + '%2.1f'%abs(fr_dur[fr_elab]).mean() + ' $\pm$ ' + ('%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_elab]))) + ' & ' + '%2.1f'%abs(fr_dur[fr_nelab]).mean() + ' $\pm$ ' + ('%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_nelab]))) + ' & ' + '%2.1f'%t4.statistic + ' & ' + ('%1.2E'%t4.pvalue)[0:4] + ' $\\times 10^{' + ('%1.2E'%t4.pvalue)[-3:] + '}$' + '\\\ \n')

        # Rope number

        f.write('Number of ropes & ' + str(len(fr_elab)) + ' & ' + str(len(fr_nelab)) + ' & . & . \\\ \n')

        f.write('\enddata' + '\n')
        f.write('\end{deluxetable*}' + '\n')

        # All good things...
        f.close()

def read(csfrm):
    '''
    Reads a flux rope footprint map for the specified data frame.
    '''

    infile = netcdf.netcdf_file(outdir + 'fr-' + csfrm + '.nc', 'r')
    frmap = infile.variables['frmap'][:,:].copy()
    #frhlcy = infile.variables['frhlcy'][:,:].copy()
    #frrext = infile.variables['frrext'][:,:].copy()
    #br0 = infile.variables['br0'][:,:].copy()
    #lat = infile.variables['lat'][:].copy()
    #lon = infile.variables['lon'][:].copy()
    infile.close()

    return frmap

if __name__ == "__main__":
    afiles = glob.glob(datdir+adatprefix+'*.nc')
    if len(afiles) == 0 : prep()
    FRoDO()
    erupt()
    plot()
    stats()
