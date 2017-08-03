# FRoDO
Flux Rope Detection and Organization

[![DOI](https://zenodo.org/badge/90994349.svg)](https://zenodo.org/badge/latestdoi/90994349)

## Goal

The FRoDO code has a few specific goals in mind:
- Automated and consistent detection of magnetic flux rope structures
- Operation independent of dataset, such that the end user only needs to provide datacubes of magnetic field components and associated grids
- Tracking of detected flux ropes in time, outputting locations and associated statistics

## Notes

For the moment this code works with a sample dataset that will be provided. This will eventually expand to include user-provided datasets, conforming with certain requirements on domain, etc. For the time being, please contact us for a sample dataset that we can provide on request.

At the moment this code works with a cadence of one day, though this will be modified and user adjustable in future releases.

## Dependencies

Python 3 is supported for the use of this code.

Several common libraries are required for this code, and can be installed (if not already) using pip:

    $ pip install matplotlib scipy numpy

A few astronomy-specific libraries prove useful:

    $ pip install astropy sunpy

## Usage

To begin the process of flux rope detection, load magnetic field netCDF data into the specified input data directory. The filename prefixes should be provided in the configuration file. Files should be labeled according to simulation day, padded with zeros to five digits. This should be followed by an underscore, and a simulation hour with two digits. For example, b_00042_12.nc .

To begin, enter into a Python 3 environment via the command line,

    $ python3

Note that for the moment, Python 3 is the recommended and supported version to execute this code. Begin by loading the FRoDO library of functions,

    >>> import FRoDO

The magnetic vector potential can then be computed for the given sample dataset using:

    >>> FRoDO.prep()

This will create a series of additional netCDF files to store the magnetic vector potential. From here, tracking can be completed with:

    >>> FRoDO.FRoDO()

The resulting output files are saved to the specified output directory. Output files are marked by 'fr-ddddd_hh.nc', and contain flux rope footprint locations, along with other associated data. Having stepped through all of the provided data, tracking is performed on detected flux rope footprints, with associated time histories stored in the subdirectory 'hist' of the specified output directory. These python Pickle files contain the lists of arrays storing the time history for each unique flux rope structure. The labeling provided in the array of footprints in the 'fr-ddddd_hh.nc' corresponds to the element of these lists containing the associated history.

Note that all large-helicity structures are tracked and stored here in this initial tracking. To filter out tall / brief features that are not quite so flux-ropish, the array fr-frg is stored in the output history directory.

To begin the process of tracking eruptive flux ropes:

    >>> FRoDO.erupt()

After processing has completed, erupting / non-erupting labels will be saved into the specified output directory, under the subdirectory 'hist'. Note that these labels must be cross-referenced with the list of flux ropes contained within fr-rfrg from earlier.

From here, the resulting time histories, eruption labels, and frame-by-frame detailed data can be recalled using the script FRoDO-read.py.

    $ python3 FRoDO-read.py

This will read in all of the aforementioned time histories, as well as eruption and radial extent/duration labels. This is mainly intended for diagnostic purposes, or to utilize as a branching off point for more customized visualizations.

To run through a more standardized set of plotting routines,

    >>> FRoDO.plot()

This will read output data, filter flux ropes accordingly, and create a standardized set of output plots for visualization. Feel free to use this as a starting point to further explore this data.

To compute some statistics for the detected flux ropes,

    >>> FRoDO.stats()

Mean values and associated standard deviations will be calculated and displayed, along with several t-tests and correlation coefficients.

Note that while all of these subroutines can be executed individually from within an interactive Python environment, FRoDO can be executed in a more script-like fashion from the command line,

    $ python3 FRoDO.py

This will run in sequence the FRoDO(), erupt(), plot(), and stats() subroutines, outputting the resulting data accordingly.
