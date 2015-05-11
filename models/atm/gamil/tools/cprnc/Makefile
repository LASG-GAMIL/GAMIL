#-----------------------------------------------------------------------
# This Makefile is for building cprnc
# on a IBM, Compaq, Linux (with PGI compiler), SGI or SUN platform
#
#
# To use set the environment variables:
#
# LIB_NETCDF --- Library directory location of NetCDF. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of NetCDF. (defaults to /usr/local/include)
#
# You also can set the environment variables:
#
# EXEDIR ------- Directory to build executable in. (Defaults to users $HOME/bin directory)
# VPATH -------- GNU make path. (Defaults to current directory)
# USER_FFLAGS -- Additional compiler flags that you want to add.
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

EXENAME := cprnc
RM      := rm

UNAMES  := $(shell uname -s)
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
FFLAGS  := -I$(INC_NETCDF) -DSUN -c -O0 -stackvar -xpp=cpp -free -g
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
  FFLAGS += -g -Mbounds -Ktrap=fp
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

# Add user defined compiler flags if set.
FFLAGS  += $(USER_FFLAGS)

shr_kind_mod.o: shr_kind_mod.F
header.o: header.F shr_kind_mod.o
cpr.o: cpr.F header.o stats.o shr_kind_mod.o nldat.o
cprtps.o: cprtps.F header.o nldat.o
header.o: header.F shr_kind_mod.o
initstats.o: initstats.F stats.o shr_kind_mod.o
ismax.o: ismax.F
lenchr.o: lenchr.F
nldat.o: nldat.F
plevs0.o: plevs0.F
prhddiff.o: prhddiff.F shr_kind_mod.o header.o
printstats.o: printstats.F stats.o shr_kind_mod.o
stats.o: stats.F shr_kind_mod.o
wrap_nf.o: shr_kind_mod.o wrap_nf.F
