# makefile

FC = gfortran
FCFLAGS = -ffree-form -ffree-line-length-none -ftree-vectorize -Wall

# FCFLAGS = -g  -O0 -ffree-line-length-none -fcheck=all -fno-check-array-temporaries -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -Wall -pedantic

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

PASTESOILCODE_OBJS = pastesoilcode.o

PASTESOIL_OBJS = pastesoil.o

NCPASTE_OBJS = ncpaste.o

PASTECOORDS_OBJS = pastecoords.o

MAKETHICKNESS_OBJS = makethickness.o

MASKLANDMASK_OBJS = masklandmask.o

MASKLANDMASK_BYTE_OBJS = masklandmask_byte.o

SOILCALC_OBJS = parametersmod.o     \
                pedotransfermod.o  \
                soilpropertiesmod.o \
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

all::	pastesoilcode ncpaste pastesoil pastecoords makethickness masklandmask masklandmask_byte soilcalc

pastesoilcode: $(PASTESOILCODE_OBJS)
	$(FC) $(FCFLAGS) -o pastesoilcode $(PASTESOILCODE_OBJS) $(LDFLAGS) $(LIBS)

pastesoil: $(PASTESOIL_OBJS)
	$(FC) $(FCFLAGS) -o pastesoil $(PASTESOIL_OBJS) $(LDFLAGS) $(LIBS)

ncpaste: $(NCPASTE_OBJS)
	$(FC) $(FCFLAGS) -o ncpaste $(NCPASTE_OBJS) $(LDFLAGS) $(LIBS)

pastecoords: $(PASTECOORDS_OBJS)
	$(FC) $(FCFLAGS) -o pastecoords $(PASTECOORDS_OBJS) $(LDFLAGS) $(LIBS)

makethickness: $(MAKETHICKNESS_OBJS)
	$(FC) $(FCFLAGS) -o makethickness $(MAKETHICKNESS_OBJS) $(LDFLAGS) $(LIBS)

masklandmask: $(MASKLANDMASK_OBJS)
	$(FC) $(FCFLAGS) -o masklandmask $(MASKLANDMASK_OBJS) $(LDFLAGS) $(LIBS)

masklandmask_byte: $(MASKLANDMASK_BYTE_OBJS)
	$(FC) $(FCFLAGS) -o masklandmask_byte $(MASKLANDMASK_BYTE_OBJS) $(LDFLAGS) $(LIBS)

soilcalc: $(SOILCALC_OBJS)
	$(FC) $(FCFLAGS) -o soilcalc $(SOILCALC_OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm *.o *.mod pastesoilcode ncpaste pastesoil pastecoords makethickness masklandmask masklandmask_byte soilcalc
