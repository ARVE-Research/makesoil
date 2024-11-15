# makefile

FC = gfortran
FCFLAGS = -ffree-form -ffree-line-length-none -ftree-vectorize -Wall

# use the command "nf-config --all" to find the location of your netCDF installation
# and enter the path next to " --prefix    ->" on the line below

netcdf = /usr/local

# should not need to modify anything below this line

# ---------------------------------------------

NC_LIB = $(netcdf)/lib
NC_INC = $(netcdf)/include

CPPFLAGS = -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)
LIBS     = -lnetcdff

# ---------------------------------------------

PASTEWRB_OBJS = pastewrb.o

PASTESOIL_OBJS = pastesoil.o

NCPASTE_OBJS = ncpaste.o

PASTECOORDS_OBJS = pastecoords.o

SOILCALC_OBJS = pedotransfermod.o \
                simplesoilmod.o   \
                soilcalc.o

# ---------------------------------------------

.SUFFIXES: .o .f90 .F90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.F90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	pastewrb ncpaste pastesoil pastecoords soilcalc

pastewrb: $(PASTEWRB_OBJS)
	$(FC) $(FCFLAGS) -o pastewrb $(PASTEWRB_OBJS) $(LDFLAGS) $(LIBS)

pastesoil: $(PASTESOIL_OBJS)
	$(FC) $(FCFLAGS) -o pastesoil $(PASTESOIL_OBJS) $(LDFLAGS) $(LIBS)

ncpaste: $(NCPASTE_OBJS)
	$(FC) $(FCFLAGS) -o ncpaste $(NCPASTE_OBJS) $(LDFLAGS) $(LIBS)

pastecoords: $(PASTECOORDS_OBJS)
	$(FC) $(FCFLAGS) -o pastecoords $(PASTECOORDS_OBJS) $(LDFLAGS) $(LIBS)

soilcalc: $(SOILCALC_OBJS)
	$(FC) $(FCFLAGS) -o soilcalc $(SOILCALC_OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm *.o *.mod pastewrb ncpaste pastesoil pastecoords soilcalc
