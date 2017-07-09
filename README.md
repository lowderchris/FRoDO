# FRoDO
Flux Rope Detection and Organization

## Goal
The goal of this project is the consolidation of various routines I've developed into a single coherent algorithm for public use and development. With that in mind, there are a few specific goals:
- Automated and consistent detection of magnetic flux rope structures
- Operation independent of dataset, such that the end user only needs to provide datacubes of magnetic field components and associated grids
- Tracking of detected flux ropes in time, outputting locations and associated statistics

## Notes

For the moment this code works with a sample dataset that will be provided. This will eventually expand to include user-provided datasets, conforming with certain requirements on domain, etc.

At the moment this code works with a cadence of one day, though this will be modified and user adjustable in future releases.

Python3 is supported for this code.

## Usage

To begin the process of flux rope detection, load netCDF data into the specified input data directory . The filename prefix should be provided in the configuration file. Files should be labeled according to simulation day, padded with zeros to five digits. From here, tracking can be completed with:

    $ python3 FRoDO.py

The resulting output files are saved to the specified output directory. Output files are marked by 'fr-#####.nc', and contain flux rope footprint locations, along with other associated data. Having stepped through all of the provided data, tracking is performed on detected flux rope footprints, with associated time histories stored in the subdirectory 'hist' of the specified output directory. These python Pickle files contain the lists of arrays storing the time history for each unique flux rope structure. The labeling provided in the array of footprints in the 'fr-#####.nc' corresponds to the element of these lists containing the associated history.

To begin the process of tracking eruptive flux ropes:

    $ python3 FRoDO-erupt.py

Note that the eruption detection code is still in the process of being ported over, and is not functional just yet...
After processing has completed, erupting / non-erupting labels will be saved into the specified output directory, under the subdirectory 'hist'.

## To-do
- [X] Create an initial README
- [ ] Create dataverse(?) repository for sample dataset
- [X] Outline of overall algorithm structure
- [ ] Break routines into modules
- [ ] Complete merging of all routines
- [X] Read parameter configuration file
- [ ] Modify for the use of differing time scales
