# FRoDO_plot.py
# Set of routines to plot more advanced plotting outputs
'''
This set of modules works alongside generated data from the FRoDO code.
Specify plotting parameters within the config.cfg file.
For all additional details, consult the aptly named README file.
'''

## To run the plot3d routine headlessly in the background,
## xvfb-run --server-args="-screen 0 1024x768x24" python3 -c "import FRoDO-plot; FRoDO-plot()"

# Import libraries
import matplotlib
matplotlib.use('Agg')

import os
import glob
import shutil
import re
import datetime
import matplotlib
from matplotlib import cm, colors
from matplotlib import pyplot
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

ref_sthresh = np.double(config['thresholds']['ref_sthresh'])
ref_ethresh = np.double(config['thresholds']['ref_ethresh'])
ref_havg = np.double(config['thresholds']['ref_havg'])

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

# Note that the plot3d routine is disabled for now, while sorting out installation issues with mayavi
def plot3d():
    '''
    Uses specified plotting parameters to animate a three-dimension visualization of magnetic flux rope evolution.
    '''

    print('I\'m still under development... check back later')

    # Import 3D plotting libraries
    from mayavi import mlab

    # Remove any existing files
    if os.path.exists(frmdir) : shutil.rmtree(frmdir)

    # Create a frame output directory
    if not os.path.exists(frmdir) : os.mkdir(frmdir)

    # Generate a list of files for animating
    files = glob.glob(outdir + 'fr-'+'*.nc')
    files.sort()

    # Read radial extent and duration filtered index
    infile = open(outdir + '/hist/fr-frg.pkl', 'rb')
    fr_frg = pickle.load(infile)
    infile.close()
    fr_frg = fr_frg.astype(np.int)

    # Read original timing data
    infile = open(outdir + '/hist/h-fr-tarr.pkl', 'rb')
    tarr = pickle.load(infile)
    infile.close()

    # Define some loop counting information
    dcount = 0
    dvcount = 0
    prntend = '\r'

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

        # Identify filtered flux rope region footpoint coordinates
        fpth = np.array([])
        fpph = np.array([])
        for reg in np.unique(frmap[frmap != 0]):
            if np.in1d(reg, fr_frg):
                frwhr = np.where(frmap == reg)
                fpth = np.append(fpth, lat[frwhr[0]])
                fpph = np.append(fpph, lon[frwhr[1]])
        fpr = np.ones(len(fpth))
        fpth = np.arcsin(-1 * fpth) + pi/2

        # Trace these flux rope fieldlines for plotting
        if len(fpr) != 0:
            fflhlcy, ffl = b.fieldlines_hlcy(fpr, fpth, fpph)
        else:
            fflhlcy = np.array([])
            ffl = np.array([])

        # Apply a final helicity filter
        havg = abs(frhlcy).mean()
        sthresh = (havg / 0.26648107322354925) * ref_sthresh
        ethresh = (havg / 0.26648107322354925) * ref_ethresh
        fl = []
        flhlcy = np.array([])
        for ifl in np.arange(len(fflhlcy)):
            if abs(fflhlcy[ifl]) > ethresh:
                fl.append(ffl[ifl])
                flhlcy = np.append(flhlcy, fflhlcy[ifl])
        del ffl
        del fflhlcy

        # Setup the scene
        sz = (1024,1024)
        if pltlght:
            bgcol = (1,1,1)
            fgcol = (0,0,0)
        else:
            bgcol = (0,0,0)
            fgcol = (1,1,1)
        if dcount != 0:
                fig.scene.disable_render = True
                mlab.clf()
                mlab.close(all=True)
        fig = mlab.figure(1, bgcolor=bgcol, fgcolor=fgcol, size=sz)
        fig.scene.disable_render = False
        mlab.view(90, 70, 9, (0,0,0))
        fig.scene.light_manager.light_mode = 'vtk'

        # Create a sun
        r = 1.0
        pi = np.pi
        cos = np.cos
        sin = np.sin
        theta, phi = np.mgrid[1:-1:180j, 0:2 * pi:360j]
        x = r * np.sin(np.arcsin(theta)+pi/2) * np.cos(phi)
        y = r * np.sin(np.arcsin(theta)+pi/2) * np.sin(phi)
        z = r * np.cos(np.arcsin(theta)+pi/2)

        sol = mlab.mesh(x, y, z, scalars=np.clip(br0,-10,10), colormap='Greys')
        sol.module_manager.scalar_lut_manager.reverse_lut = True

        frnorm = matplotlib.colors.Normalize(vmin=-3,vmax=3)
        frsmap = matplotlib.cm.ScalarMappable(norm=frnorm, cmap='RdBu_r')

        cval_hlcy = np.array([abs(flhlcy.min()),abs(flhlcy.max())]).min()
        cmap_hlcy = cm.get_cmap('RdBu_r')
        norm_hlcy = matplotlib.colors.Normalize(vmin=-cval_hlcy, vmax=cval_hlcy)
        smap_hlcy = matplotlib.cm.ScalarMappable(norm=norm_hlcy, cmap=cmap_hlcy)

        # Trace a set of background fieldlines
        nph = 12
        nth = 12
        r = 2.5
        fth1 = np.linspace(0.05*np.pi, 0.95*np.pi, nth)
        fr1 = np.array([r])
        fph1 = np.linspace(0, 2*np.pi, nph, endpoint=False)
        fr0, fth0, fph0 = np.meshgrid(fr1, fth1, fph1)
        fr0 = fr0.flatten()
        fth0 = fth0.flatten()
        fph0 = fph0.flatten()
        ahlcy, afl = b.fieldlines_hlcy(fr0, fth0, fph0)
        for f in np.arange(len(afl)):
                if np.mod(f,50) == 0:
                    print('Plotting background fieldline ' + '%05.f'%f + '/' + '%05.f'%len(afl))
                slen = np.sqrt((afl[f][0] * sin(afl[f][1]) * cos(afl[f][2]))**2 +            (afl[f][0] * sin(afl[f][1]) * sin(afl[f][2]))**2 + (afl[f][0] * cos(afl[f][1]))**2)
                tfl = mlab.plot3d(afl[f][0] * sin(afl[f][1]) * cos(afl[f][2]), afl[f][0] *   sin(afl[f][1]) * sin(afl[f][2]), afl[f][0] * cos(afl[f][1]), slen, tube_radius=None,         colormap='YlGn', vmin=1.0, vmax=2.5, reset_zoom=False, opacity=flalpha)
                tfl.module_manager.scalar_lut_manager.reverse_lut = True

        fig.scene.disable_render = True

        # Trace flux rope fieldlines
        fr_count = np.int64(0)
        for i in np.arange(len(fl)):
                flcol = frsmap.to_rgba(flhlcy[i])
                if np.mod(fr_count,500) == 0:
                    print('Plotting flux rope fieldline ' + '%05.f'%fr_count + '/' + '%05.f'%len(fl))
                mlab.plot3d(fl[i][0] * sin(fl[i][1]) * cos(fl[i][2]), fl[i][0] *     sin(fl[i][1]) * sin(fl[i][2]), fl[i][0] * cos(fl[i][1]), tube_radius=None, color = flcol[0:3], reset_zoom=False, opacity=flalpha)
                fr_count = fr_count + 1

        fig.scene.disable_render = False
        if pltanno:
            mlab.text(0.01,0.01,csfrm + ' ' + tarr[dcount].strftime('%Y-%m-%d %H:%M'),opacity=0.5,width=0.07)

        # Save output frame(s) and render additional frames
        for iv in np.arange(1):
            mlab.view(vph[dvcount],vth[dvcount],vr[dvcount], (0,0,0))
            mlab.savefig(datdir + 'scratch/'+frmdir+'frm'+'%05.f'%dvcount+'.png')
            dvcount = dvcount + 1

        # Diagnostic readouts!
        time1 = datetime.datetime.now()
        if dcount == 0:
            timedel = (time1 - time0)
        else:
            timedel = ((time1 - time0) + timedel) / 2
        timeeta = (nfrm - (dcount+1)) * timedel + time1
        if dcount == (nfrm - 1) : prntend = '\n'
        print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta), end=prntend)

        dcount = dcount + 1

    # Animate
    animate('plot3d')

def plot2d():
    '''
    Flattens data to a two-dimensional representation for a latitude-longitude projection plot of flux rope fieldlines.
    '''

    # A quick announcement
    print('Plotting detected flux ropes in two dimensions...')

    # Remove any existing files
    if os.path.exists(frmdir) : shutil.rmtree(frmdir)

    # Create a frame output directory
    if not os.path.exists(frmdir) : os.mkdir(frmdir)

    # Create an output plot directory if needed
    if not os.path.exists('plt') : os.mkdir('plt')

    # Generate a list of files for animating
    files = glob.glob(outdir + 'fr-'+'*.nc')
    files.sort()
    nfrm = len(files)

    # Read radial extent and duration filtered index
    infile = open(outdir + '/hist/fr-frg.pkl', 'rb')
    fr_frg = pickle.load(infile)
    infile.close()
    fr_frg = fr_frg.astype(np.int)

    # Read original timing data
    infile = open(outdir + '/hist/h-fr-tarr.pkl', 'rb')
    tarr = pickle.load(infile)
    infile.close()

    # Define some loop counting information
    dcount = 0
    prntend = '\r'

    # Iterate through these files
    for file in files:

        # Define some timing
        time0 = datetime.datetime.now()
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

        # Identify filtered flux rope region footpoint coordinates
        fpth = np.array([])
        fpph = np.array([])
        for reg in np.unique(frmap[frmap != 0]):
            if np.in1d(reg, fr_frg):
                frwhr = np.where(frmap == reg)
                fpth = np.append(fpth, lat[frwhr[0]])
                fpph = np.append(fpph, lon[frwhr[1]])
        fpr = np.ones(len(fpth))
        fpth = np.arcsin(-1 * fpth) + pi/2

        # Trace these flux rope fieldlines for plotting
        if len(fpr) != 0:
            fflhlcy, ffl = b.fieldlines_hlcy(fpr, fpth, fpph)
        else:
            fflhlcy = np.array([])
            ffl = np.array([])

        # Apply a final helicity filter
        havg = abs(frhlcy).mean()
        sthresh = (havg / 0.26648107322354925) * ref_sthresh
        ethresh = (havg / 0.26648107322354925) * ref_ethresh
        fl = []
        flhlcy = np.array([])
        for ifl in np.arange(len(fflhlcy)):
            if abs(fflhlcy[ifl]) > ethresh:
                fl.append(ffl[ifl])
                flhlcy = np.append(flhlcy, fflhlcy[ifl])
        del ffl
        del fflhlcy

        # Begin plotting
        # Setup color mapping
        if dcount == 0:
            frnorm = matplotlib.colors.Normalize(vmin=-3,vmax=3)
            frsmap = matplotlib.cm.ScalarMappable(norm=frnorm, cmap='RdBu_r')
            lnsty = 'solid'
            lnthick = 0.5
            lnalpha = 0.5

        # Plot everything
        pyplot.close()
        f, (ax) = pyplot.subplots(1, figsize=[7,3.5])
        ax.contour(lon, lat, br0, levels=[-20,-10,-5,5,10,20], colors='k', linewidths=lnthick)
        for i in np.arange(len(fl)):
            flcol = frsmap.to_rgba(flhlcy[i])

            alfpr1 = fl[i][0][0]
            alfpr2 = fl[i][0][-1]
            alfpth1 = np.sin(fl[i][1][0] + pi/2)
            alfpth2 = np.sin(fl[i][1][-1] + pi/2)
            alfpph1 = fl[i][2][0]
            alfpph2 = fl[i][2][-1]

            arr = np.roll(fl[i][2],1) - fl[i][2]
            cph = (np.where(abs(arr) > 6.0))[0]

            spt = 0
            ept = len(fl[i][2])-1
            if np.size(cph) != 0:
                for dpt in cph:
                    ax.plot(fl[i][2][spt:dpt], np.sin(fl[i][1][spt:dpt]+pi/2), color=flcol, ls=lnsty, alpha=lnalpha, lw=lnthick)
                    spt = dpt
                ax.plot(fl[i][2][dpt:ept], np.sin(fl[i][1][dpt:ept]+pi/2), color=flcol, ls=lnsty, alpha=lnalpha, lw=lnthick)
            else:
                ax.plot(fl[i][2], np.sin(fl[i][1]+pi/2), color=flcol, ls=lnsty, alpha=lnalpha, lw=lnthick)

        ax.set_title(csfrm.replace('_', '\_') + ' ' + tarr[dcount].strftime('%Y-%m-%d %H:%M'))
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Sine latitude')
        ax.set_xlim([0,2*pi])
        ax.set_ylim([-1,1])
        pyplot.tight_layout()
        pyplot.savefig(frmdir + 'frm'+'%05.f'%dcount+'.png')

        # Diagnostic readouts!
        time1 = datetime.datetime.now()
        if dcount == 0:
            timedel = (time1 - time0)
        else:
            timedel = ((time1 - time0) + timedel) / 2
        timeeta = (nfrm - (dcount+1)) * timedel + time1
        if dcount == (nfrm - 1) : prntend = '\n'
        print('Frame ' + '%05.f'%(dcount+1) + ' / ' + '%05.f'%nfrm + ' - ' + str(timedel) + 's - ETA ' + str(timeeta), end=prntend)

        dcount = dcount + 1

    # Animate
    animate('plot2d')

def animate(filename, frmrt=10):
   os.system('ffmpeg -r ' + str(frmrt) + ' -i ' + frmdir + 'frm%05d.png -vcodec libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -q 0 ./plt/' + filename + '.mp4') 
