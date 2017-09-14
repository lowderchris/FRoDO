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

    $ pip install numpy scipy matplotlib

A few astronomy-specific libraries prove useful:

    $ pip install astropy sunpy

A library for image processing

    $ pip install scikit-image

A library for better colormap management:

    $ pip install palettable

A library for reading configuration files:

    $ pip install configparser

See footnotes for a full list of dependency versions that have been tested and confirmed.

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

Please note that there may be issues with compatibility between Python 2 and Python 3 usage with certain versions of f2py.

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

The resulting output files are saved to the specified output directory. Output files are marked with the prefix 'fr-', and contain flux rope footprint locations, along with other associated data. Having stepped through all of the provided data, tracking is performed on detected flux rope footprints, with associated time histories stored in the subdirectory 'hist' of the specified output directory. Flux rope quantities are stored within the fr-hist.nc file, with full histories stored within .pkl files. The labeling provided in the array of footprints corresponds to the element of these lists containing the associated history.

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

## FRoDO_plot

The two and three dimensional plotting routines outlined in this section are still under development, and require package dependencies that are... tricky to install and update. Proceed with caution all ye who enter here...

### Dependencies

The Mayavi set of routines is the primary dependency for these three dimensional plotting routines, and also the source of most of the headaches involved. [Installation instructions](http://docs.enthought.com/mayavi/mayavi/installation.html) are available, where using a recommended Python bundle (Anaconda, Enthought, etc) is recommended. The author has had luck in the past with installation on a Macintosh machine with the methodology,

    brew install qt
    brew install vtk --with-python3 --without-python
    pip3 install mayavi

Note that to enable off-screen plotting, the xvfb package is also required. This should be available on some Macintosh machines, and should also be available through the [XQuartz](https://www.xquartz.org) tools.

The two dimensional plotting routine is much simpler in scope, and relies only on standard matplotlib plotting libraries.

### Usage

Before starting, a few specified parameters must be set under the plot3d section of config.cfg. These parameters specify a temporary output frame directory, viewing angles, and other plotting toggles. To begin with plotting, enter a Python environment and import the plot3d libraries,

    $ python3
    >>> import FRoDO_plot

To execute the two-dimensional animation process,

    >>> FRoDO_plot.plot2d()

This will animate frames individually, outputting to the specified output frame directory. These frames will then be animated in an output file plt/plot2d.mp4

To run through the three-dimensional animation process, execute this with the command,

    >>> FRoDO_plot.plot3d()

(Note that the moment that this routine is most easily compatible with the Python 2 version of the Mayavi libraries.). This will render individual frames, outputting data to the specified output frame directory. On completion, these frames will be animated into the file plt/plot3d.mp4

After generating a series of frames with these two routines, an animation can be created manually using,

    >>> FRoDO_plot.animate(filename, frmrt=30)

Here the filename can be specified, along with the output integer framerate.

## Footnotes

Note that for the above packages, the following versions have been tested and confirmed:

- astropy (2.0.1)
- cycler (0.10.0)
- decorator (4.1.2)
- matplotlib (2.0.2)
- networkx (1.11)
- numpy (1.13.1)
- olefile (0.44)
- palettable (3.0.0)
- pandas (0.20.3)
- Pillow (4.2.1)
- pip (9.0.1)
- py (1.4.34)
- pyparsing (2.2.0)
- pytest (3.2.1)
- python-dateutil (2.6.1)
- pytz (2017.2)
- PyWavelets (0.5.2)
- scikit-image (0.13.0)
- scipy (0.19.1)
- setuptools (36.2.7)
- six (1.10.0)
- sunpy (0.8.0)
- wheel (0.29.0)
