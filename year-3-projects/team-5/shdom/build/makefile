# Makefile for SHDOM and related programs in distribution
#
# Set the UseMPI and UseNetcdf variables to make the desired shdom90 version.

FC     = gfortran     # your Fortran 90 compiler here (same as used by mpif90)
FFLAGS = -O3     #   your flags here
MPIFC  = mpif90  # MPI Fortran 90 compiling command for shdom90
MPIFLAGS = -O3   #   MPI compiling flags
LIBS   = 

SHDOMOBJS   = shdomsub1.o shdomsub2.o shdomsub3.o fftpack.o ocean_brdf.o 
SHDOM90OBJS = $(SHDOMOBJS)

UseMPI = no
UseNetcdf = no

NetcdfHome = /usr/local/netcdf-3.6.3
NetcdfInl  = -I$(NetcdfHome)/include
NetcdfLib  = -L$(NetcdfHome)/lib -lnetcdf


ifeq ($(UseMPI),no)
  SHDOM90OBJS += shdom_nompi.o
  F90 = $(FC)
  F90FLAGS = $(FFLAGS)
else
  SHDOM90OBJS += shdom_mpi.o
  F90 = $(MPIFC)
  F90FLAGS = $(MPIFLAGS)
endif

ifeq ($(UseNetcdf),no)
  SHDOM90OBJS += shdom_nonetcdf.o
else
  SHDOM90OBJS += shdom_netcdf.o
  F90FLAGS += $(NetcdfInl)
  LIBS   += $(NetcdfLib)
endif

all: shdom90 shdom propgen make_mie_table make_ice_table plotscattab \
      cloudprp ckdswrrtm ckdlwrrtm ckdfu 

shdom: shdom.o shdom_nompi.o $(SHDOMOBJS)
	$(FC) $(FFLAGS) -o shdom shdom.o shdom_nompi.o $(SHDOMOBJS) 

shdom90: shdom90.o $(SHDOM90OBJS)
	$(F90) $(F90FLAGS) -o shdom90 shdom90.o $(SHDOM90OBJS) $(LIBS)

shdom90.o: shdom90.f90 $(SHDOM90OBJS)
	$(F90) -c $(F90FLAGS) shdom90.f90

shdom_mpi.o: shdom_mpi.f90
	$(F90) -c $(F90FLAGS) shdom_mpi.f90

shdom_nonetcdf.o: shdom_nonetcdf.f90
	$(F90) -c $(F90FLAGS) shdom_nonetcdf.f90

shdom_netcdf.o: shdom_netcdf.f90
	$(F90) -c $(F90FLAGS) shdom_netcdf.f90


AsciiPropertyFileToNetcdf: AsciiPropertyFileToNetcdf.o 
	$(F90) $(F90FLAGS) -o AsciiPropertyFileToNetcdf AsciiPropertyFileToNetcdf.o shdom_netcdf.o  $(LIBS)
AsciiPropertyFileToNetcdf.o: AsciiPropertyFileToNetcdf.f90 shdom_netcdf.o
	$(F90) -c $(F90FLAGS) AsciiPropertyFileToNetcdf.f90

propgen: propgen.f90
	$(FC) $(FFLAGS) propgen.f90  -o propgen

make_ice_table: make_ice_table.f90
	$(FC) $(FFLAGS) make_ice_table.f90  -o make_ice_table

make_mie_table: make_mie_table.f90
	$(FC) $(FFLAGS) make_mie_table.f90 mieindsub.f  -o make_mie_table

plotscattab: plotscattab.f90
	$(FC) $(FFLAGS) plotscattab.f90  -o plotscattab


cloudprp: cloudprp.o mieindsub.o
	$(FC) $(FFLAGS) cloudprp.o mieindsub.o  -o cloudprp

ckdswrrtm: ckdswrrtm.o swrrtm_kdata.o
	$(FC) $(FFLAGS) ckdswrrtm.o swrrtm_kdata.o  -o ckdswrrtm

ckdlwrrtm: ckdlwrrtm.o lwrrtm_kdata.o
	$(FC) $(FFLAGS) ckdlwrrtm.o lwrrtm_kdata.o  -o ckdlwrrtm

ckdfu: ckdfu.o fukdata.o
	$(FC) $(FFLAGS) ckdfu.o fukdata.o  -o ckdfu


.f.o :; $(FC) -c $(FFLAGS) $*.f

