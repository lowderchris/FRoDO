# FRoDO-stats.py
# Script to calculate a few statistics on the resulting flux rope data

# Import libraries

import scipy
from scipy import stats
import pickle

import configparser

# Read parameters from a configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

datdir = config['paths']['datdir']
datprefix = config['paths']['datprefix']
outdir = config['paths']['outdir']

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

# Calculate and print some statistics

print('Mean erupting unsigned net helicity (Mx^2)', '%1.2E'%abs(fr_nhlcy[fr_elab]).mean())
print('Mean non-erupting unsigned net helicity (Mx^2)', '%1.2E'%abs(fr_nhlcy[fr_nelab]).mean())
print('Standard deviation of erupting unsigned net helicity (Mx^2)', '%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_elab])))
print('Standard deviation non-erupting unsigned net helicity (Mx^2)', '%1.2E'%scipy.stats.tstd(abs(fr_nhlcy[fr_nelab])))

print('Mean erupting unsigned magnetic flux (Mx)','%1.2E'%abs(fr_uflux[fr_elab]).mean())
print('Mean non-erupting unsigned magnetic flux (Mx)','%1.2E'%abs(fr_uflux[fr_nelab]).mean())
print('Standard deviation erupting unsigned magnetic flux (Mx)','%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_elab])))
print('Standard deviation non-erupting unsigned magnetic flux (Mx)','%1.2E'%scipy.stats.tstd(abs(fr_uflux[fr_nelab])))

print('Mean erupting footprint area (cm^2)','%1.2E'%abs(fr_area[fr_elab]).mean())
print('Mean non-erupting footprint area (cm^2)','%1.2E'%abs(fr_area[fr_nelab]).mean())
print('Standard deviation erupting footprint area (cm^2)','%1.2E'%scipy.stats.tstd(abs(fr_area[fr_elab])))
print('Standard deviation non-erupting footprint area (cm^2)','%1.2E'%scipy.stats.tstd(abs(fr_area[fr_nelab])))

print('Mean erupting duration (days)','%2.1f'%abs(fr_dur[fr_elab]).mean())
print('Mean non-erupting duration (days)','%2.1f'%abs(fr_dur[fr_nelab]).mean())
print('Standard deviation erupting duration (days)','%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_elab])))
print('Standard deviation non-erupting duration (days)','%2.1f'%scipy.stats.tstd(abs(fr_dur[fr_nelab])))

# Calculate the t-test for the erupting and non-erupting sets
# The independent t-test should be used here, as these are separate sets
# Note that Welchs t-test is used here, as we have differing sample sizes, with unsure equality of population variance.
# This calculation returns the calculated t-statistic value, and the two-tailed p-value.

print('T-test for unsigned net helicity:', scipy.stats.ttest_ind(abs(fr_nhlcy[fr_elab]), abs(fr_nhlcy[fr_nelab]),equal_var=False))

print('T-test for unsigned magnetic flux:', scipy.stats.ttest_ind(fr_uflux[fr_elab], fr_uflux[fr_nelab],equal_var=False))

print('T-test for footprint area:', scipy.stats.ttest_ind(fr_area[fr_elab], fr_area[fr_nelab],equal_var=False))

print('T-test for duration:', scipy.stats.ttest_ind(fr_dur[fr_elab], fr_dur[fr_nelab],equal_var=False))

# Calculate linear fits for scatter distributions, along with appropriate statistics.

# Duration and unsigned net helicity
print('Spearman rank-order calculation for non-erupting duration and unsigned net helicity:',scipy.stats.spearmanr(fr_dur[fr_nelab], abs(fr_nhlcy[fr_nelab])))
print('Pearson correlation coefficient for non-erupting duration and unsigned net helicity:', scipy.stats.pearsonr(fr_dur[fr_nelab], abs(fr_nhlcy[fr_nelab])))

print('Spearman rank-order calculation for erupting duration and unsigned net helicity:',scipy.stats.spearmanr(fr_dur[fr_elab], abs(fr_nhlcy[fr_elab])))
print('Pearson correlation coefficient for erupting duration and unsigned net helicity:', scipy.stats.pearsonr(fr_dur[fr_elab], abs(fr_nhlcy[fr_elab])))

# Duration and unsigned magnetic flux
print('Spearman rank-order calculation for non-erupting duration and unsigned magnetic flux:',scipy.stats.spearmanr(fr_dur[fr_nelab], fr_uflux[fr_nelab]))
print('Pearson correlation coefficient for non-erupting duration and unsigned magnetic flux:', scipy.stats.pearsonr(fr_dur[fr_nelab], fr_uflux[fr_nelab]))

print('Spearman rank-order calculation for erupting duration and unsigned magnetic flux:',scipy.stats.spearmanr(fr_dur[fr_elab], fr_uflux[fr_elab]))
print('Pearson correlation coefficient for erupting duration and unsigned magnetic flux:', scipy.stats.pearsonr(fr_dur[fr_elab], fr_uflux[fr_elab]))

# Unsigned magnetic flux and unsigned net helicity
print('Spearman rank-order calculation for non-erupting unsigned magnetic flux  and unsigned net helicity:',scipy.stats.spearmanr(fr_uflux[fr_nelab], abs(fr_nhlcy[fr_nelab])))
print('Pearson correlation coefficient for non-erupting unsigned magnetic flux and unsigned net helicity:', scipy.stats.pearsonr(fr_uflux[fr_nelab], abs(fr_nhlcy[fr_nelab])))

print('Spearman rank-order calculation for erupting unsigned magnetic flux  and unsigned net helicity:',scipy.stats.spearmanr(fr_uflux[fr_elab], abs(fr_nhlcy[fr_elab])))
print('Pearson correlation coefficient for erupting unsigned magnetic flux and unsigned net helicity:', scipy.stats.pearsonr(fr_uflux[fr_elab], abs(fr_nhlcy[fr_elab])))
