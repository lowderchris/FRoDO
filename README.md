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

## Dependencies

Python 3 is supported for the use of this code.

Several common libraries are required for this code, and can be installed (if not already) using pip:

    $ pip install matplotlib scipy numpy

A few astronomy-specific libraries prove useful:

    $ pip install astropy sunpy

### Compiled dependencies

Two sets of FORTRAN codes must first be compiled for use with the main FRoDO Python code.

A fieldline tracer FORTRAN code, tracer.f90, can be compiled directly with the commands,

    $ cd tracer
    $ f2py -m tracer --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c tracer.f90

A set of routines for the calculation of the magnetic vector potential can be computed using the appropriate commands,

    $ cd compA
    $ f2py -m compA --fcompiler=gfortran -I/path/to/gfortran/modules/ -L/path/to/libnetcdff.so.5 -lnetcdff --f90flags='-fopenmp' -lgomp -c main.f90 shared.f90 input_output.f90 vectorpot.f90

Note that this requires specification of the location of netCDF libraries on the current machine. This may vary for the particular machine and method of installation of netCDF libraries. As an example, for a Linux machine, libraries are located using,

    -I/usr/lib64/gfortran/modules -L/usr/lib64/libnetcdff.so.5 -lnetcdff

For a sample Macintosh machine (with Homebrew netCDF), these should be located with:

    -I/usr/local/include/ -L/usr/lib/libnetcdff.so.5 -lnetcdff

## Data

For the current version of FRoDO, input magnetic field data files must be in the netCDF format, with filenames in chronological order, and have the following variables contained within:

- (r, th, ph) - Coordinate arrays
- (br, bth, bph) - Datacubes of magnetic field spherical components
- date - Datestring in the format %Y%b%d_%H%M

## Usage

To begin the process of flux rope detection, load magnetic field netCDF data into the specified input data directory. The filename prefixes should be provided in the configuration file. 

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

To run through a standardized set of plotting routines,

    >>> FRoDO.plot()

This will read output data, filter flux ropes accordingly, and create a standardized set of output plots for visualization. Feel free to use this as a starting point to further explore this data.

To compute some statistics for the detected flux ropes,

    >>> FRoDO.stats()

Mean values and associated standard deviations will be calculated, along with several t-tests and correlation coefficients. This data will be written to the file stats.txt, located in the user-specified output directory. An optional flag, FRoDO.stats(tex=True) will in addition output stats-tex.tex, two pre-built deluxetable snippets that display relevant flux rope statistics.

Note that while all of these subroutines can be executed individually from within an interactive Python environment, FRoDO can be executed in a more script-like fashion from the command line,

    $ python3 FRoDO.py

This will run in sequence the prep(), FRoDO(), erupt(), plot(), and stats() subroutines, outputting the resulting data accordingly. Note that the prep() function will only be executed if the magnetic vector potential has not already been calculated, with appropriate data in the input directory.
