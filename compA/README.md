# compA

To compile the existing FORTRAN routines to be used internally within the FRoDO scripts,

For a Linux machine:
f2py -m compA --fcompiler=gfortran -I/usr/lib64/gfortran/modules -L/usr/lib64/libnetcdff.so.5 -lnetcdff --f90flags='-fopenmp' -lgomp -c main.f90 shared.f90 input_output.f90 vectorpot.f90

For a Macintosh machine:
f2py -m compA --fcompiler=gfortran -I/usr/local/include/ -L/usr/lib/libnetcdff.so.5 -lnetcdff --f90flags='-fopenmp' -lgomp -c main.f90 shared.f90 input_output.f90 vectorpot.f90

Note that this requires specification of the location of netCDF libraries on the current machine. For one of the linux machines, these paths should be:
-I/usr/lib64/gfortran/modules -L/usr/lib64/libnetcdff.so.5 -lnetcdff
For a macintosh machine, these should be located with:
-I/usr/local/include/ -L/usr/lib/libnetcdff.so.5 -lnetcdff

This will compile the FORTRAN routines for internal use within python.
