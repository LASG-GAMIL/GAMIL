#-----------------------------------------------------------------------
# This Makefile is for building cprlndnc
# on a CRAY, SGI or SUN platform
#
#------------------------------------------------------------------------
# Set up special characters
null  :=


#
# Default variables
#

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

EXENAME := cprlndnc
RM      := rm

UNAMES  := $(shell uname -s)
UNAMEM  := $(findstring CRAY,$(shell uname -m))
# If path to source code not given
ifeq ($(VPATH),$(null))
  VPATH:= .
endif
dirs := $(subst :,$(space),$(VPATH))

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.c $(dir)/*.F $(dir)/*.F90)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

# If executable directory not given
ifeq ($(EXEDIR),$(null))
  EXEDIR  := $(HOME)/bin
  ifeq ($(UNAMES),IRIX64)
    EXEDIR  := $(HOME)/bin/IRIX64
  endif
endif

.SUFFIXES:
.SUFFIXES: .F .o

$(EXEDIR)/$(EXENAME):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------
ifeq ($(UNAMEM),CRAY)
FC      := f90

CPP      := cpp
FPPFLAGS := -DCRAY -N -P
CFLAGS   := -DCRAY -g
FFLAGS   := -DCRAY -c -G1 -Otask0 -I$(INC_NETCDF)
LDFLAGS  := -L $(HOME)/lib,$(LIB_NETCDF) -lncaro,timing,netcdf
MODULES   =

.F.o:
	$(CPP) $(FPPFLAGS) $< > $*.f90
	$(FC) $(FFLAGS) $*.f90
#	$(RM) $*.f90

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)
FC      := f90

CPP     := /lib/cpp
FFLAGS  := -64 -DSGI -I$(INC_NETCDF) -mips4 -g -c -C -freeform
LDFLAGS := -64 -g -L$(LIB_NETCDF) -lnetcdf -lfpe
MODULES = *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------
ifeq ($(UNAMES),OSF1)
FC      := f90

FFLAGS  := -I$(INC_NETCDF) -g -c -free
LDFLAGS := -L$(LIB_NETCDF) -lnetcdf
MODULES = *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#-----------------------------------------------------------------------
# SUN
#-----------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC      := f90
CPP     := cpp
CFLAGS  := -DSUN -I.. -g
FFLAGS  := -I$(INC_NETCDF) -DSUN -c -O2 -stackvar -xpp=cpp -free
LDFLAGS := -L$(LIB_NETCDF) -L$(HOME)/lib -lnetcdf
MODULES = *.M *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# LINUX
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC    := pgf90

FFLAGS := -DLINUX -I$(INC_NETCDF) -c -Mfreeform -Mrecursive -Mdalign -Mfreeform
ifeq ($(DEBUG),TRUE)
  FFLAGS += -g -Mbounds -Ktrap=fp -DDEBUG
else
  FFLAGS += -fast 
endif

LDFLAGS := -g -L$(LIB_NETCDF) -lnetcdf

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC    := xlf90

FFLAGS  := -I$(INC_NETCDF) -c -qarch=auto -qstrict -O3 -qdpc=e
LDFLAGS := -bmaxdata:0x20000000 -L$(LIB_NETCDF) -lnetcdf


.F.o:
	$(FC) $(FFLAGS) $<

endif

#
# Clean rule
#
clean:
	$(RM) -f $(OBJS) $(MODULES) $(EXEDIR)/$(EXENAME)

shr_kind_mod.o: shr_kind_mod.F
header.o: header.F shr_kind_mod.o
cpr.o: cpr.F header.o stats.o shr_kind_mod.o nldat.o
cprtps.o: cprtps.F header.o nldat.o
header.o: header.F shr_kind_mod.o
initstats.o: initstats.F stats.o shr_kind_mod.o
ismax.o: ismax.F
lenchr.o: lenchr.F
nldat.o: nldat.F
prhddiff.o: prhddiff.F shr_kind_mod.o header.o
printstats.o: printstats.F stats.o shr_kind_mod.o
stats.o: stats.F shr_kind_mod.o
wrap_nf.o: shr_kind_mod.o wrap_nf.F
