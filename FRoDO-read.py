# FRoDO-read.py
# Script to read processed data from FRoDO

import pickle
from scipy.io import netcdf

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

datdir = config['paths']['datdir']
outdir = config['paths']['outdir']

sfrm = np.int(config['times']['sfrm'])
efrm = np.int(config['times']['efrm'])
dfrm = np.int(config['times']['dfrm'])

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

# Read a particular frame for reference
cfrm = 899
csfrm = '%05.f'%cfrm

infile = netcdf.netcdf_file(outdir + 'fr-' + csfrm + '.nc', 'r')
frmap = infile.variables['frmap'][:,:].copy()
frhlcy = infile.variables['frhlcy'][:,:].copy()
frrext = infile.variables['frrext'][:,:].copy()
br0 = infile.variables['br0'][:,:].copy()
lat = infile.variables['lat'][:].copy()
lon = infile.variables['lon'][:].copy()
infile.close()
