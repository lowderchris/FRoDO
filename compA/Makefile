# Set FORTRAN90 compiler
FC = gfortran  # for Beltrami

# Location of files for netcdf library
NETCDF = -I /usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64/libnetcdff.so.5 -lnetcdff

# Set compiler flags for optimized gfortran
FFLAGS = -O3 -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5 -fopenmp  #-fbounds-check

# Set build parameters
TARGET = computeA

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin

OBJFILES = shared.o input_output.o vectorpot.o main.o
FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

	
# All the dependencies
shared.o: shared.f90
vectorpot.o: vectorpot.f90 shared.o
input_output.o: input_output.f90 shared.o
main.o: main.f90 input_output.o vectorpot.o shared.o
