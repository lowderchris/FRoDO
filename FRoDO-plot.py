# FRoDO-plot.py
# Script to visualize detected flux rope data

# Import libraries
import os
import numpy as np
import scipy
import pickle
from scipy.io import netcdf

from matplotlib import gridspec
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
import mpl_toolkits.axisartist.floating_axes as floating_axes

import configparser

# Define color tables
import palettable
cols = (palettable.colorbrewer.get_map('Paired', 'Qualitative', 12)).mpl_colors

# Define some plotting parameters
fscale = 1.0    # Relative figure scale
fnapp = ''      # Label to append to plot filename
defcol = 'k'    # Color
legfsz = 8      # Legend font size
gsize = 1e20    # Butterfly glyph scaling size

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

datdir = config['paths']['datdir']
outdir = config['paths']['outdir']

sfrm = np.int(config['times']['sfrm'])
efrm = np.int(config['times']['efrm'])
dfrm = np.int(config['times']['dfrm'])

# Create output directories if needed
os.system("mkdir " + 'plt')

# Remove any existing files
os.system("rm " + 'plt/' + "*")

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
fr_etarr = fr_etarr.astype(np.int)

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

# Compute the flux and helicity ejection rates
ejr_uflux = np.zeros(len(tarr), dtype=np.double)
ejr_nhlcy = np.zeros(len(tarr), dtype=np.double)
for i in fr_efpt:
    ejr_uflux[fr_time[i].astype(np.int)] = ejr_uflux[fr_time[i].astype(np.int)] + fr_uflux[i]
    ejr_nhlcy[fr_time[i].astype(np.int)] = ejr_nhlcy[fr_time[i].astype(np.int)] + abs(fr_nhlcy[i])

# Begin plotting

# Radial graph of eruption mean latitudes
fpmlatne = np.histogram(np.arcsin(fr_mlat[fr_nelab])*180./pi, bins=90, range=[-90,90], density=False)
fpmlate = np.histogram(np.arcsin(fr_mlat[fr_elab])*180./pi, bins=90, range=[-90,90], density=False)
whg = (180.) / len(fpmlatne[0])
rscl = np.array([fpmlatne[0].max(), fpmlate[0].max()]).max() * 1.2

f = figure(figsize=fscale*np.array([3.35,2.45]))
ax, aax = setup_polaraxis(f, 111, rscl)
aax.plot(fpmlatne[1][0:-1], fpmlatne[0], label='Non-erupting')
aax.plot(fpmlate[1][0:-1], fpmlate[0], label='Erupting')
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
sctne = ax1.scatter(netarr, np.arcsin(fr_mlat[fr_nelab])*180./pi, c=fr_nhlcy[fr_nelab], cmap='RdBu_r', s=fr_area[fr_nelab]/gsize,edgecolors='None',vmin=-0.5e43,vmax=0.5e43, alpha=1.0)
scte = ax2.scatter(etarr, np.arcsin(fr_mlat[fr_elab])*180./pi, c=fr_nhlcy[fr_elab], cmap='RdBu_r', s=fr_area[fr_elab]/gsize,edgecolors='None',vmin=-0.5e43,vmax=0.5e43, alpha=1.0)
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
f.savefig('plt/fr-bfly-nshlcy-erupt'+fnapp+'.pdf')

# A few hexbin scatter plots

# Scatter of duration and net helicity magnitude
ne_coef = np.polyfit(np.log10(fr_dur[fr_nelab]), np.log10(abs(fr_nhlcy[fr_nelab])),1)
ne_polyn = np.poly1d(ne_coef)
ne_xs = np.linspace(0.1,2.8)
ne_ys = ne_polyn(ne_xs)

e_coef = np.polyfit(np.log10(fr_dur[fr_elab]), np.log10(abs(fr_nhlcy[fr_elab])),1)
e_polyn = np.poly1d(e_coef)
e_xs = np.linspace(0.1,2.8)
e_ys = e_polyn(e_xs)

f, (ax1,ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
hb1 = ax1.hexbin(fr_dur[fr_nelab], abs(fr_nhlcy[fr_nelab]), xscale='log', yscale='log', gridsize=25, extent=[0,3,40,44], cmap='Blues', linewidths=0)
hb2 = ax2.hexbin(fr_dur[fr_elab], abs(fr_nhlcy[fr_elab]), xscale='log', yscale='log', gridsize=25, extent=[0,3,40,44], cmap='Oranges', linewidths=0)
ax2.set_xlabel('Duration [days]')
ax2.set_ylabel('$|$H$|$ [Mx$^2$]')
ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+'} x^{'+'%.03f'%ne_coef[0]+'}$')
ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+'} x^{'+'%.03f'%e_coef[0]+'}$')
leg1 = ax1.legend(loc=4, fontsize=legfsz)
leg2 = ax2.legend(loc=4, fontsize=legfsz)
ax1.set_xlim(1e0,1e3)
ax1.set_ylim(1e40,1e44)
ax2.set_xlim(1e0,1e3)
ax2.set_ylim(1e40,1e44)
ax1.xaxis.set_ticklabels([])
colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
tight_layout()
savefig('plt/fr-sct-dur-nhlcy'+fnapp+'.pdf')

# Scatter of duration and unsigned magnetic flux
ne_coef = np.polyfit(np.log10(fr_dur[fr_nelab]), np.log10(fr_uflux[fr_nelab]),1)
ne_polyn = np.poly1d(ne_coef)
ne_xs = np.linspace(0.1,2.8)
ne_ys = ne_polyn(ne_xs)

e_coef = np.polyfit(np.log10(fr_dur[fr_elab]), np.log10(fr_uflux[fr_elab]),1)
e_polyn = np.poly1d(e_coef)
e_xs = np.linspace(0.1,2.8)
e_ys = e_polyn(e_xs)

f, (ax1,ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
hb1 = ax1.hexbin(fr_dur[fr_nelab], fr_uflux[fr_nelab], xscale='log', yscale='log', gridsize=25, extent=[0,3,18,23], cmap='Blues', linewidths=0)
hb2 = ax2.hexbin(fr_dur[fr_elab], fr_uflux[fr_elab], xscale='log', yscale='log', gridsize=25, extent=[0,3,18,23], cmap='Oranges', linewidths=0)
ax2.set_xlabel('Duration [days]')
ax2.set_ylabel('$|\Phi_m|$ [Mx]')
ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+  '} x^{'+'%.03f'%ne_coef[0]+'}$')
ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+  '} x^{'+'%.03f'%e_coef[0]+'}$')
leg1 = ax1.legend(loc=4, fontsize=legfsz)
leg2 = ax2.legend(loc=4, fontsize=legfsz)
ax1.set_xlim(1e0,1e3)
ax1.set_ylim(1e18,1e23)
ax2.set_xlim(1e0,1e3)
ax2.set_ylim(1e18,1e23)
ax1.xaxis.set_ticklabels([])
colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
tight_layout()
savefig('plt/fr-sct-dur-uflux'+fnapp+'.pdf')

# Scatter of unsigned magnetic flux and net helicity magnitude
ne_coef = np.polyfit(np.log10(fr_uflux[fr_nelab]), np.log10(abs(fr_nhlcy[fr_nelab])),1)
ne_polyn = np.poly1d(ne_coef)
ne_xs = np.linspace(18.1,22.8)
ne_ys = ne_polyn(ne_xs)

e_coef = np.polyfit(np.log10(fr_uflux[fr_elab]), np.log10(abs(fr_nhlcy[fr_elab])),1)
e_polyn = np.poly1d(e_coef)
e_xs = np.linspace(18.1,22.8)
e_ys = e_polyn(e_xs)

f, (ax1, ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
hb1=ax1.hexbin(fr_uflux[fr_nelab], abs(fr_nhlcy[fr_nelab]), xscale='log', yscale='log', gridsize=25, cmap='Blues', extent=[18,23,40,44], linewidths=0)
hb2=ax2.hexbin(fr_uflux[fr_elab], abs(fr_nhlcy[fr_elab]), xscale='log', yscale='log', gridsize=25, cmap='Oranges', extent=[18,23,40,44], alpha=1.0, linewidths=0)
#ax1.grid()
#ax2.grid()
ax2.set_xlabel(r'$|\Phi_m|$ [Mx]')
ax2.set_ylabel(r'$|$H$|$ [Mx$^2$]')
ax1.plot(10**ne_xs, 10**ne_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%ne_coef[1]+  '} x^{'+'%.03f'%ne_coef[0]+'}$')
ax2.plot(10**e_xs, 10**e_ys, color='0.5', linestyle='dashed', label='$10^{'+'%.01f'%e_coef[1]+  '} x^{'+'%.03f'%e_coef[0]+'}$')
leg1 = ax1.legend(loc=4, fontsize=legfsz)
leg2 = ax2.legend(loc=4, fontsize=legfsz)
ax1.set_xlim(1e18,1e23)
ax1.set_ylim(1e40,1e44)
ax2.set_xlim(1e18,1e23)
ax2.set_ylim(1e40,1e44)
ax1.xaxis.set_ticklabels([])
colorbar(hb1,ax=ax1, label='N [bin$^{-1}$]')
colorbar(hb2,ax=ax2, label='N [bin$^{-1}$]')
tight_layout()
f.savefig('plt/fr-sct-uflux-nhlcy'+fnapp+'.pdf')

# Helicity and unsigned flux ejection rates
f, (ax1, ax2) = subplots(2,1, figsize=fscale*np.array([3.35,5]))
ax1.plot(tarr, ejr_uflux, label='1-day')
ax1.plot(tarr, scipy.convolve(ejr_uflux, np.ones(7)/7., mode='same'), label='7-day', color=defcol)
#ax1.plot(tarr, scipy.convolve(ejr_uflux, np.ones(27)/27., mode='same'), label='27-day', color=defcol)
#ax1.plot(tarr, scipy.convolve(ejr_uflux, np.ones(164)/164., mode='same'), label='6-month', color=cols[1])
ax2.plot(tarr, ejr_nhlcy, label='1-day')
ax2.plot(tarr, scipy.convolve(ejr_nhlcy, np.ones(7)/7., mode='same'), label='7-day', color=defcol)
#ax2.plot(tarr, scipy.convolve(ejr_nhlcy, np.ones(27)/27., mode='same'), label='27-day', color=defcol)
#ax2.plot(tarr, scipy.convolve(ejr_nhlcy, np.ones(164)/164., mode='same'), label='6-month', color=cols[1])

ax2.set_xlabel('Date')
ax1.set_ylabel(r'$\Phi_m$ / day [Mx day$^{-1}$]')
ax2.set_ylabel(r'H / day [Mx$^2$ day$^{-1}$]')
ax1.xaxis.set_ticklabels([])
ax1.legend()
for label in ax2.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
f.autofmt_xdate() # Optional toggle to sort out overlapping dates
tight_layout()
savefig('plt/fr_ejr'+fnapp+'.pdf')
